// src/emmer_generator.cpp
#include "emmer_generator.hpp"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

//=============================================================================
// EmmerGeneratorConfig constructor
//=============================================================================

EmmerGeneratorConfig::EmmerGeneratorConfig() {
    model_index = 0;
    align_output = 1;
    set_refmac_blur = 1;
    blur = 0.0;
    symmetry_expansion = 1;
    cutoff_level = 1e-5;
    rate = 1.5;
}

//=============================================================================
// EmmerGaussianCoeff methods
//=============================================================================

double EmmerGaussianCoeff::fe0() const {
    double sum = 0.0;
    for (int i = 0; i < EMMER_N_GAUSS; i++) {
    sum += a[i];
    }
    return sum;
}

//=============================================================================
// EmmerFFT methods
//=============================================================================

void EmmerFFT::init() {
    nx = ny = nz = 0;
    nreals = ncomplex = 0;
    real_in = nullptr;
    complex_out = nullptr;
    kernel = nullptr;
    fft_desc = nullptr;
}


void EmmerFFT::cleanup() {
    if (real_in) {
    ippsFree(real_in);
    real_in = nullptr;
    }
    if (complex_out) {
    ippsFree(complex_out);
    complex_out = nullptr;
    }
    if (kernel) {
    ippsFree(kernel);
    kernel = nullptr;
    }
    if (fft_desc) {
    DftiFreeDescriptor(&fft_desc);
    fft_desc = nullptr;
    }
    nx = ny = nz = 0;
    nreals = ncomplex = 0;
}

int EmmerFFT::setup(int nx_, int ny_, int nz_) {
    nx = nx_;
    ny = ny_;
    nz = nz_;
    nreals = (int64_t)nx * ny * nz;
    ncomplex = (int64_t)nz * ny * (nx / 2 + 1);

    // Allocate buffers
    real_in = (Ipp64f*)ippsMalloc_64f(nreals);
    complex_out = (Ipp64fc*)ippsMalloc_64fc(ncomplex);
    kernel = (Ipp64fc*)ippsMalloc_64fc(ncomplex);

    if (!real_in || !complex_out || !kernel) {
    cleanup();
    return -1;
    }

    // Create FFT descriptor
    MKL_LONG dims[3] = { nz, ny, nx };
    MKL_LONG status = DftiCreateDescriptor(&fft_desc, DFTI_DOUBLE, DFTI_REAL, 3, dims);
    if (status != 0) return -2;

    status = DftiSetValue(fft_desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    if (status != 0) return -3;

    status = DftiSetValue(fft_desc, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    if (status != 0) return -4;

    // Set strides for real data (row-major order: [z][y][x])
    MKL_LONG rstrides[4] = { 0, ny * nx, nx, 1 };
    status = DftiSetValue(fft_desc, DFTI_INPUT_STRIDES, rstrides);
    if (status != 0) return -5;

    // Set strides for complex data (conjugate-even storage)
    MKL_LONG cstrides[4] = { 0, ny * (nx / 2 + 1), (nx / 2 + 1), 1 };
    status = DftiSetValue(fft_desc, DFTI_OUTPUT_STRIDES, cstrides);
    if (status != 0) return -6;

    status = DftiCommitDescriptor(fft_desc);
    if (status != 0) return -7;

    return 0;
}

int EmmerFFT::forward(Ipp64f* input) {
    if (!fft_desc || !input || !complex_out) return -1;

    // Copy input to real_in
    ippsCopy_64f(input, real_in, nreals);

    // Forward FFT (real-to-complex)
    MKL_LONG status = DftiComputeForward(fft_desc, real_in, complex_out);
    return (status == 0) ? 0 : -2;
}

int EmmerFFT::applyUnblur(double blur, double grid_spacing) {
    if (!fft_desc || !complex_out || !kernel) return -1;

    // Copy complex_out to kernel
    ippsCopy_64fc(complex_out, kernel, ncomplex);

    // Apply unblur: multiply by exp(-blur/4 * d*^2)
    // d*^2 = (h^2/a^2 + k^2/b^2 + l^2/c^2) in reciprocal space
    // For grid points, the frequency coordinates are:
    // For i from 0 to nx-1, frequency = i/(nx * grid_spacing)

    double inv_nx = 1.0 / (nx * grid_spacing);
    double inv_ny = 1.0 / (ny * grid_spacing);
    double inv_nz = 1.0 / (nz * grid_spacing);

#pragma omp parallel for collapse(2)
    for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
    // For each (k,j), i runs 0 to nx/2 (due to conjugate-even storage)
    int64_t base_idx = ((int64_t)k * ny + j) * (nx / 2 + 1);

    // Calculate frequency components
    double fz = (k <= nz / 2) ? k * inv_nz : (k - nz) * inv_nz;
    double fy = (j <= ny / 2) ? j * inv_ny : (j - ny) * inv_ny;

    for (int i = 0; i <= nx / 2; i++) {
    double fx = i * inv_nx;

    // d*^2 = fx^2 + fy^2 + fz^2
    double dst2 = fx * fx + fy * fy + fz * fz;

    // Unblur factor = exp(-blur/4 * d*^2)
    double factor = exp(-blur * 0.25 * dst2);

    int64_t idx = base_idx + i;
    kernel[idx].re *= (Ipp64f)factor;
    kernel[idx].im *= (Ipp64f)factor;
    }
    }
    }

    return 0;
}

int EmmerFFT::inverse(Ipp64f* output) {
    if (!fft_desc || !kernel || !output) return -1;

    // Inverse FFT (complex-to-real)
    MKL_LONG status = DftiComputeBackward(fft_desc, kernel, output);
    return (status == 0) ? 0 : -2;
}

//=============================================================================
// Complete coefficients from GEMMI's c4322.hpp (International Tables Vol. C)
//=============================================================================

const EmmerGaussianCoeff EmmerGenerator::COEFFICIENTS[] = {
    // Format: {symbol, {a1, a2, a3, a4, a5}, {b1, b2, b3, b4, b5}, Z}
    // H (Z=1)
    {"H",   {0.0349, 0.1201, 0.1970, 0.0573, 0.1195},
    {0.5347, 3.5867, 12.3471, 18.9525, 38.6269}, 1},
    // He (Z=2)
    {"He",  {0.0317, 0.0838, 0.1526, 0.1334, 0.0164},
    {0.2507, 1.4751, 4.4938, 12.6646, 31.1653}, 2},
    // Li (Z=3)
    {"Li",  {0.0750, 0.2249, 0.5548, 1.4954, 0.9354},
    {0.3864, 2.9383, 15.3829, 53.5545, 138.7337}, 3},
    // Be (Z=4)
    {"Be",  {0.0780, 0.2210, 0.6740, 1.3867, 0.6925},
    {0.3131, 2.2381, 10.1517, 30.9061, 78.3273}, 4},
    // B (Z=5)
    {"B",   {0.0909, 0.2551, 0.7738, 1.2136, 0.4606},
    {0.2995, 2.1155, 8.3816, 24.1292, 63.1314}, 5},
    // C (Z=6)
    {"C",   {0.0893, 0.2563, 0.7570, 1.0487, 0.3575},
    {0.2465, 1.7100, 6.4094, 18.6113, 50.2523}, 6},
    // N (Z=7)
    {"N",   {0.1022, 0.3219, 0.7982, 0.8197, 0.1715},
    {0.2451, 1.7481, 6.1925, 17.3894, 48.1431}, 7},
    // O (Z=8)
    {"O",   {0.0974, 0.2921, 0.6910, 0.6990, 0.2039},
    {0.2067, 1.3815, 4.6943, 12.7105, 32.4726}, 8},
    // F (Z=9)
    {"F",   {0.1083, 0.3175, 0.6487, 0.5846, 0.1421},
    {0.2057, 1.3439, 4.2788, 11.3932, 28.7881}, 9},
    // Ne (Z=10)
    {"Ne",  {0.1269, 0.3535, 0.5582, 0.4674, 0.1460},
    {0.2200, 1.3779, 4.0203, 9.4934, 23.1278}, 10},
    // Na (Z=11)
    {"Na",  {0.2142, 0.6853, 0.7692, 1.6589, 1.4482},
    {0.3334, 2.3446, 10.0830, 48.3037, 138.2700}, 11},
    // Mg (Z=12)
    {"Mg",  {0.2314, 0.6866, 0.9677, 2.1882, 1.1339},
    {0.3278, 2.2720, 10.9241, 39.2898, 101.9748}, 12},
    // Al (Z=13)
    {"Al",  {0.2390, 0.6573, 1.2011, 2.5586, 1.2312},
    {0.3138, 2.1063, 10.4163, 34.4552, 98.5344}, 13},
    // Si (Z=14)
    {"Si",  {0.2519, 0.6372, 1.3795, 2.5082, 1.0500},
    {0.3075, 2.0174, 9.6746, 29.3744, 80.4732}, 14},
    // P (Z=15)
    {"P",   {0.2548, 0.6106, 1.4541, 2.3204, 0.8477},
    {0.2908, 1.8740, 8.5176, 24.3434, 63.2996}, 15},
    // S (Z=16)
    {"S",   {0.2497, 0.5628, 1.3899, 2.1865, 0.7715},
    {0.2681, 1.6711, 7.0267, 19.5377, 50.3888}, 16},
    // Cl (Z=17)
    {"Cl",  {0.2443, 0.5397, 1.3919, 2.0197, 0.6621},
    {0.2468, 1.5242, 6.1537, 16.6687, 42.3086}, 17},
    // Ar (Z=18)
    {"Ar",  {0.2385, 0.5017, 1.3428, 1.8899, 0.6079},
    {0.2289, 1.3694, 5.2561, 14.0928, 35.5361}, 18},
    // K (Z=19)
    {"K",   {0.4115, 1.4031, 2.2784, 2.6742, 2.2162},
    {0.3703, 3.3874, 13.1029, 68.9592, 194.4329}, 19},
    // Ca (Z=20)
    {"Ca",  {0.4054, 1.3880, 2.1602, 3.7532, 2.2063},
    {0.3499, 3.0991, 11.9608, 53.9353, 142.3892}, 20},
    // Sc (Z=21)
    {"Sc",  {0.3787, 1.2181, 2.0594, 3.2618, 2.3870},
    {0.3133, 2.5856, 9.5813, 41.7688, 116.7282}, 21},
    // Ti (Z=22)
    {"Ti",  {0.3825, 1.2598, 2.0008, 3.0617, 2.0694},
    {0.3040, 2.4863, 9.2783, 39.0751, 109.4583}, 22},
    // V (Z=23)
    {"V",   {0.3876, 1.2750, 1.9109, 2.8314, 1.8979},
    {0.2967, 2.3780, 8.7981, 35.9528, 101.7201}, 23},
    // Cr (Z=24)
    {"Cr",  {0.4046, 1.3696, 1.8941, 2.0800, 1.2196},
    {0.2986, 2.3958, 9.1406, 37.4701, 113.7121}, 24},
    // Mn (Z=25)
    {"Mn",  {0.3796, 1.2094, 1.7815, 2.5420, 1.5937},
    {0.2699, 2.0455, 7.4726, 31.0604, 91.5622}, 25},
    // Fe (Z=26)
    {"Fe",  {0.3946, 1.2725, 1.7031, 2.3140, 1.4795},
    {0.2717, 2.0443, 7.6007, 29.9714, 86.2265}, 26},
    // Co (Z=27)
    {"Co",  {0.4118, 1.3161, 1.6493, 2.1930, 1.2830},
    {0.2742, 2.0372, 7.7205, 29.9680, 84.9383}, 27},
    // Ni (Z=28)
    {"Ni",  {0.3860, 1.1765, 1.5451, 2.0730, 1.3814},
    {0.2478, 1.7660, 6.3107, 25.2204, 74.3146}, 28},
    // Cu (Z=29)
    {"Cu",  {0.4314, 1.3208, 1.5236, 1.4671, 0.8562},
    {0.2694, 1.9223, 7.3474, 28.9892, 90.6246}, 29},
    // Zn (Z=30)
    {"Zn",  {0.4288, 1.2646, 1.4472, 1.8294, 1.0934},
    {0.2593, 1.7998, 6.7500, 25.5860, 73.5284}, 30},
    // Ga (Z=31)
    {"Ga",  {0.4818, 1.4032, 1.6561, 2.4605, 1.1054},
    {0.2825, 1.9785, 8.7546, 32.5238, 98.5523}, 31},
    // Ge (Z=32)
    {"Ge",  {0.4655, 1.3014, 1.6088, 2.6998, 1.3003},
    {0.2647, 1.7926, 7.6071, 26.5541, 77.5238}, 32},
    // As (Z=33)
    {"As",  {0.4517, 1.2229, 1.5852, 2.7958, 1.2638},
    {0.2493, 1.6436, 6.8154, 22.3681, 62.0390}, 33},
    // Se (Z=34)
    {"Se",  {0.4477, 1.1678, 1.5843, 2.8087, 1.1956},
    {0.2405, 1.5442, 6.3231, 19.4610, 52.0233}, 34},
    // Br (Z=35)
    {"Br",  {0.4798, 1.1948, 1.8695, 2.6953, 0.8203},
    {0.2504, 1.5963, 6.9653, 19.8492, 50.3233}, 35},
    // Kr (Z=36)
    {"Kr",  {0.4546, 1.0993, 1.7696, 2.7068, 0.8672},
    {0.2309, 1.4279, 5.9449, 16.6752, 42.2243}, 36},
    // Rb (Z=37)
    {"Rb",  {1.0160, 2.8528, 3.5466, -7.7804, 12.1148},
    {0.4853, 5.0925, 25.7851, 130.4515, 138.6775}, 37},
    // Sr (Z=38)
    {"Sr",  {0.6703, 1.4926, 3.3368, 4.4600, 3.1501},
    {0.3190, 2.2287, 10.3504, 52.3291, 151.2216}, 38},
    // Y (Z=39)
    {"Y",   {0.6894, 1.5474, 3.2450, 4.2126, 2.9764},
    {0.3189, 2.2904, 10.0062, 44.0771, 125.0120}, 39},
    // Zr (Z=40)
    {"Zr",  {0.6719, 1.4684, 3.1668, 3.9557, 2.8920},
    {0.3036, 2.1249, 8.9236, 36.8458, 108.2049}, 40},
    // Nb (Z=41)
    {"Nb",  {0.6123, 1.2677, 3.0348, 3.3841, 2.3683},
    {0.2709, 1.7683, 7.2489, 27.9465, 98.5624}, 41},
    // Mo (Z=42)
    {"Mo",  {0.6773, 1.4798, 3.1788, 3.0824, 1.8384},
    {0.2920, 2.0606, 8.1129, 30.5336, 100.0658}, 42},
    // Tc (Z=43)
    {"Tc",  {0.7082, 1.6392, 3.1993, 3.4327, 1.8711},
    {0.2976, 2.2106, 8.5246, 33.1456, 96.6377}, 43},
    // Ru (Z=44)
    {"Ru",  {0.6735, 1.4934, 3.0966, 2.7254, 1.5597},
    {0.2773, 1.9716, 7.3249, 26.6891, 90.5581}, 44},
    // Rh (Z=45)
    {"Rh",  {0.6413, 1.3690, 2.9854, 2.6952, 1.5433},
    {0.2580, 1.7721, 6.3854, 23.2549, 85.1517}, 45},
    // Pd (Z=46)
    {"Pd",  {0.5904, 1.1775, 2.6519, 2.2875, 0.8689},
    {0.2324, 1.5019, 5.1591, 15.5428, 46.8213}, 46},
    // Ag (Z=47)
    {"Ag",  {0.6377, 1.3790, 2.8294, 2.3631, 1.4553},
    {0.2466, 1.6974, 5.7656, 20.0943, 76.7372}, 47},
    // Cd (Z=48)
    {"Cd",  {0.6364, 1.4247, 2.7802, 2.5973, 1.7886},
    {0.2407, 1.6823, 5.6588, 20.7219, 69.1109}, 48},
    // In (Z=49)
    {"In",  {0.6768, 1.6589, 2.7740, 3.1835, 2.1326},
    {0.2522, 1.8545, 6.2936, 25.1457, 84.5448}, 49},
    // Sn (Z=50)
    {"Sn",  {0.7224, 1.9610, 2.7161, 3.5603, 1.8972},
    {0.2651, 2.0604, 7.3011, 27.5493, 81.3349}, 50},
    // Sb (Z=51)
    {"Sb",  {0.7106, 1.9247, 2.6149, 3.8322, 1.8899},
    {0.2562, 1.9646, 6.8852, 24.7648, 68.9168}, 51},
    // Te (Z=52)
    {"Te",  {0.6947, 1.8690, 2.5356, 4.0013, 1.8955},
    {0.2459, 1.8542, 6.4411, 22.1730, 59.2206}, 52},
    // I (Z=53)
    {"I",   {0.7047, 1.9484, 2.5940, 4.1526, 1.5057},
    {0.2455, 1.8638, 6.7639, 21.8007, 56.4395}, 53},
    // Xe (Z=54)
    {"Xe",  {0.6737, 1.7908, 2.4129, 4.2100, 1.7058},
    {0.2305, 1.6890, 5.8218, 18.3928, 47.2496}, 54},
    // Cs (Z=55)
    {"Cs",  {1.2704, 3.8018, 5.6618, 0.9205, 4.8105},
    {0.4356, 4.2058, 23.4342, 136.7783, 171.7561}, 55},
    // Ba (Z=56)
    {"Ba",  {0.9049, 2.6076, 4.8498, 5.1603, 4.7388},
    {0.3066, 2.4363, 12.1821, 54.6135, 161.9978}, 56},
    // La (Z=57)
    {"La",  {0.8405, 2.3863, 4.6139, 5.1514, 4.7949},
    {0.2791, 2.1410, 10.3400, 41.9148, 132.0204}, 57},
    // Ce (Z=58)
    {"Ce",  {0.8551, 2.3915, 4.5772, 5.0278, 4.5118},
    {0.2805, 2.1200, 10.1808, 42.0633, 130.9893}, 58},
    // Pr (Z=59)
    {"Pr",  {0.9096, 2.5313, 4.5266, 4.6376, 4.3690},
    {0.2939, 2.2471, 10.8266, 48.8842, 147.6020}, 59},
    // Nd (Z=60)
    {"Nd",  {0.8807, 2.4183, 4.4448, 4.6858, 4.1725},
    {0.2802, 2.0836, 10.0357, 47.4506, 146.9976}, 60},
    // Pm (Z=61)
    {"Pm",  {0.9471, 2.5463, 4.3523, 4.4789, 3.9080},
    {0.2977, 2.2276, 10.5762, 49.3619, 145.3580}, 61},
    // Sm (Z=62)
    {"Sm",  {0.9699, 2.5837, 4.2778, 4.4575, 3.5985},
    {0.3003, 2.2447, 10.6487, 50.7994, 146.4179}, 62},
    // Eu (Z=63)
    {"Eu",  {0.8694, 2.2413, 3.9196, 3.9694, 4.5498},
    {0.2653, 1.8590, 8.3998, 36.7397, 125.7089}, 63},
    // Gd (Z=64)
    {"Gd",  {0.9673, 2.4702, 4.1148, 4.4972, 3.2099},
    {0.2909, 2.1014, 9.7067, 43.4270, 125.9474}, 64},
    // Tb (Z=65)
    {"Tb",  {0.9325, 2.3673, 3.8791, 3.9674, 3.7996},
    {0.2761, 1.9511, 8.9296, 41.5937, 131.0122}, 65},
    // Dy (Z=66)
    {"Dy",  {0.9505, 2.3705, 3.8218, 4.0471, 3.4451},
    {0.2773, 1.9469, 8.8862, 43.0938, 133.1396}, 66},
    // Ho (Z=67)
    {"Ho",  {0.9248, 2.2428, 3.6182, 3.7910, 3.7912},
    {0.2660, 1.8183, 7.9655, 33.1129, 101.8139}, 67},
    // Er (Z=68)
    {"Er",  {1.0373, 2.4824, 3.6558, 3.8925, 3.0056},
    {0.2944, 2.0797, 9.4156, 45.8056, 132.7720}, 68},
    // Tm (Z=69)
    {"Tm",  {1.0075, 2.3787, 3.5440, 3.6932, 3.1759},
    {0.2816, 1.9486, 8.7162, 41.8420, 125.0320}, 69},
    // Yb (Z=70)
    {"Yb",  {1.0347, 2.3911, 3.4619, 3.6556, 3.0052},
    {0.2855, 1.9679, 8.7619, 42.3304, 125.6499}, 70},
    // Lu (Z=71)
    {"Lu",  {0.9927, 2.2436, 3.3554, 3.7813, 3.0994},
    {0.2701, 1.8073, 7.8112, 34.4849, 103.3526}, 71},
    // Hf (Z=72)
    {"Hf",  {1.0295, 2.2911, 3.4110, 3.9497, 2.4925},
    {0.2761, 1.8625, 8.0961, 34.2712, 98.5295}, 72},
    // Ta (Z=73)
    {"Ta",  {1.0190, 2.2291, 3.4097, 3.9252, 2.2679},
    {0.2694, 1.7962, 7.6944, 31.0942, 91.1089}, 73},
    // W (Z=74)
    {"W",   {0.9853, 2.1167, 3.3570, 3.7981, 2.2798},
    {0.2569, 1.6745, 7.0098, 26.9234, 81.3910}, 74},
    // Re (Z=75)
    {"Re",  {0.9914, 2.0858, 3.4531, 3.8812, 1.8526},
    {0.2548, 1.6518, 6.8845, 26.7234, 81.7215}, 75},
    // Os (Z=76)
    {"Os",  {0.9813, 2.0322, 3.3665, 3.6235, 1.9741},
    {0.2487, 1.5973, 6.4737, 23.2817, 70.9254}, 76},
    // Ir (Z=77)
    {"Ir",  {1.0194, 2.0645, 3.4425, 3.4914, 1.6976},
    {0.2554, 1.6475, 6.5966, 23.2269, 70.0272}, 77},
    // Pt (Z=78)
    {"Pt",  {0.9148, 1.8096, 3.2134, 3.2953, 1.5754},
    {0.2263, 1.3813, 5.3243, 17.5987, 60.0171}, 78},
    // Au (Z=79)
    {"Au",  {0.9674, 1.8916, 3.3993, 3.0524, 1.2607},
    {0.2358, 1.4712, 5.6758, 18.7119, 61.5286}, 79},
    // Hg (Z=80)
    {"Hg",  {1.0033, 1.9469, 3.4396, 3.1548, 1.4180},
    {0.2413, 1.5298, 5.8009, 19.4520, 60.5753}, 80},
    // Tl (Z=81)
    {"Tl",  {1.0689, 2.1038, 3.6039, 3.4927, 1.8283},
    {0.2540, 1.6715, 6.3509, 23.1531, 78.7099}, 81},
    // Pb (Z=82)
    {"Pb",  {1.0891, 2.1867, 3.6160, 3.8031, 1.8994},
    {0.2552, 1.7174, 6.5131, 23.9170, 74.7039}, 82},
    // Bi (Z=83)
    {"Bi",  {1.1007, 2.2306, 3.5689, 4.1549, 2.0382},
    {0.2546, 1.7351, 6.4948, 23.6464, 70.3780}, 83},
    // Po (Z=84)
    {"Po",  {1.1568, 2.4353, 3.6459, 4.4064, 1.7179},
    {0.2648, 1.8786, 7.1749, 25.1766, 69.2821}, 84},
    // At (Z=85)
    {"At",  {1.0909, 2.1976, 3.3831, 4.6700, 2.1277},
    {0.2466, 1.6707, 6.0197, 20.7657, 57.2663}, 85},
    // Rn (Z=86)
    {"Rn",  {1.0756, 2.1630, 3.3178, 4.8852, 2.0489},
    {0.2402, 1.6169, 5.7644, 19.4568, 52.5009}, 86},
    // Fr (Z=87)
    {"Fr",  {1.4282, 3.5081, 5.6767, 4.1964, 3.8946},
    {0.3183, 2.6889, 13.4816, 54.3866, 200.8321}, 87},
    // Ra (Z=88)
    {"Ra",  {1.3127, 3.1243, 5.2988, 5.3891, 5.4133},
    {0.2887, 2.2897, 10.8276, 43.5389, 145.6109}, 88},
    // Ac (Z=89)
    {"Ac",  {1.3128, 3.1021, 5.3385, 5.9611, 4.7562},
    {0.2861, 2.2509, 10.5287, 41.7796, 128.2973}, 89},
    // Th (Z=90)
    {"Th",  {1.2553, 2.9178, 5.0862, 6.1206, 4.7122},
    {0.2701, 2.0636, 9.3051, 34.5977, 107.9200}, 90},
    // Pa (Z=91)
    {"Pa",  {1.3218, 3.1444, 5.4371, 5.6444, 4.0107},
    {0.2827, 2.2250, 10.2454, 41.1162, 124.4449}, 91},
    // U (Z=92)
    {"U",   {1.3382, 3.2043, 5.4558, 5.4839, 3.6342},
    {0.2838, 2.2452, 10.2519, 41.7251, 124.9023}, 92},
    // Np (Z=93)
    {"Np",  {1.5193, 4.0053, 6.5327, -0.1402, 6.7489},
    {0.3213, 2.8206, 14.8878, 68.9103, 81.7257}, 93},
    // Pu (Z=94)
    {"Pu",  {1.3517, 3.2937, 5.3213, 4.6466, 3.5714},
    {0.2813, 2.2418, 9.9952, 42.7939, 132.1739}, 94},
    // Am (Z=95)
    {"Am",  {1.2135, 2.7962, 4.7545, 4.5731, 4.4786},
    {0.2483, 1.8437, 7.5421, 29.3841, 112.4579}, 95},
    // Cm (Z=96)
    {"Cm",  {1.2937, 3.1100, 5.0393, 4.7546, 3.5031},
    {0.2638, 2.0341, 8.7101, 35.2992, 109.4972}, 96},
    // Bk (Z=97)
    {"Bk",  {1.2915, 3.1023, 4.9309, 4.6009, 3.4661},
    {0.2611, 2.0023, 8.4377, 34.1559, 105.8911}, 97},
    // Cf (Z=98)
    {"Cf",  {1.2089, 2.7391, 4.3482, 4.0047, 4.6497},
    {0.2421, 1.7487, 6.7262, 23.2153, 80.3108}, 98}
};

const int EmmerGenerator::N_COEFFICIENTS = sizeof(COEFFICIENTS) / sizeof(EmmerGaussianCoeff);

//=============================================================================
// EmmerGenerator constructor/destructor
//=============================================================================

EmmerGenerator::EmmerGenerator() {
    m_atoms = nullptr;
    m_nAtoms = 0;
    m_output = nullptr;
    m_work_grid = nullptr;
    m_nElements = 0;
    m_nx = m_ny = m_nz = 0;
    m_origin[0] = m_origin[1] = m_origin[2] = 0.0;
    m_grid_spacing = 0.0;
    m_d_min = 0.0;
    m_apix[0] = m_apix[1] = m_apix[2] = 0.0;

    // Initialize element data
    for (int i = 0; i < EMMER_MAX_ELEMENTS; i++) {
        m_element_data[i].element_idx = -1;
        m_element_data[i].atom_count = 0;
        m_element_data[i].atom_indices = nullptr;
        m_element_data[i].avg_b_factor = 0.0;
        for (int g = 0; g < EMMER_N_GAUSS; g++) {
            m_element_data[i].precalc[g].amplitude = 0.0;
            m_element_data[i].precalc[g].width = 0.0;
            m_element_data[i].precalc[g].radius = 0.0;
        }
    }

    m_fft.init();  // Make sure this initializes all FFT members to nullptr
}

EmmerGenerator::~EmmerGenerator() {
    cleanup();
}

void EmmerGenerator::cleanup() {
    // Free atom index arrays
    for (int i = 0; i < EMMER_MAX_ELEMENTS; i++) {
    if (m_element_data[i].atom_indices) {
    free(m_element_data[i].atom_indices);
    m_element_data[i].atom_indices = nullptr;
    }
    // Reset other fields to prevent use-after-free
    m_element_data[i].element_idx = -1;
    m_element_data[i].atom_count = 0;
    m_element_data[i].avg_b_factor = 0.0;
    }
    m_nElements = 0;

    // Free working grid with null check
    if (m_work_grid) {
    ippsFree(m_work_grid);
    m_work_grid = nullptr;
    }

    // Free output map with null check
    if (m_output) {
    ippsFree(m_output);
    m_output = nullptr;
    }

    // Clean up FFT (which has its own null checks)
    m_fft.cleanup();

    // Reset scalar members to prevent use after destruction
    m_atoms = nullptr;
    m_nAtoms = 0;
    m_nx = m_ny = m_nz = 0;
    m_origin[0] = m_origin[1] = m_origin[2] = 0.0;
    m_grid_spacing = 0.0;
    m_d_min = 0.0;
    m_apix[0] = m_apix[1] = m_apix[2] = 0.0;
}

//=============================================================================
// Private methods
//=============================================================================

int EmmerGenerator::findUniqueElements() {
    if (!m_atoms || m_nAtoms <= 0) return -1;

    // Temporary arrays to count elements
    int element_counts[EMMER_MAX_ELEMENTS] = { 0 };
    int element_present[EMMER_MAX_ELEMENTS] = { 0 };

    // First pass: count atoms per element
    for (int i = 0; i < m_nAtoms; i++) {
    int Z = AtomUtils::get_atomic_number(m_atoms[i].element);
    if (Z <= 0 || Z > EMMER_MAX_ELEMENTS) {
    Z = 6;  // Default to carbon (Z=6)
    }
    element_counts[Z - 1]++;
    element_present[Z - 1] = 1;
    }

    // Count unique elements
    m_nElements = 0;
    for (int z = 0; z < EMMER_MAX_ELEMENTS; z++) {
    if (element_present[z]) {
    m_nElements++;
    }
    }

    if (m_nElements == 0) return -2;

    // Allocate atom index arrays and fill element data
    int elem_idx = 0;
    for (int z = 0; z < EMMER_MAX_ELEMENTS; z++) {
    if (element_present[z]) {
    m_element_data[elem_idx].element_idx = z;  // 0-based index in COEFFICIENTS
    m_element_data[elem_idx].atom_count = element_counts[z];
    m_element_data[elem_idx].atom_indices = (int*)malloc(element_counts[z] * sizeof(int));
    if (!m_element_data[elem_idx].atom_indices) {
    // Clean up already allocated arrays
    for (int j = 0; j < elem_idx; j++) {
    if (m_element_data[j].atom_indices) {
    free(m_element_data[j].atom_indices);
    m_element_data[j].atom_indices = nullptr;
    }
    }
    return -3;
    }
    elem_idx++;
    }
    }

    // Second pass: fill atom indices
    int* current_pos = (int*)malloc(m_nElements * sizeof(int));
    if (!current_pos) {
    return -4;
    }

    for (int e = 0; e < m_nElements; e++) {
    current_pos[e] = 0;
    }

    for (int i = 0; i < m_nAtoms; i++) {
    int Z = AtomUtils::get_atomic_number(m_atoms[i].element);
    if (Z <= 0 || Z > EMMER_MAX_ELEMENTS) {
    Z = 6;
    }

    // Find element data for this Z
    for (int e = 0; e < m_nElements; e++) {
    if (m_element_data[e].element_idx == Z - 1) {
    int pos = current_pos[e];
    if (pos < m_element_data[e].atom_count) {
    m_element_data[e].atom_indices[pos] = i;
    current_pos[e]++;
    }
    break;
    }
    }
    }

    free(current_pos);
    return m_nElements;
}

int EmmerGenerator::resolveParameters(const Ipp64f* cell_dimensions) {
    if (cell_dimensions) {
    // We have unit cell dimensions
    m_nx = (int)ceil(cell_dimensions[0] / m_grid_spacing);
    m_ny = (int)ceil(cell_dimensions[1] / m_grid_spacing);
    m_nz = (int)ceil(cell_dimensions[2] / m_grid_spacing);

    m_origin[0] = 0.0;
    m_origin[1] = 0.0;
    m_origin[2] = 0.0;
    }
    else {
    // No unit cell - calculate bounding box from atoms
    Ipp64f xmin = m_atoms[0].x, xmax = m_atoms[0].x;
    Ipp64f ymin = m_atoms[0].y, ymax = m_atoms[0].y;
    Ipp64f zmin = m_atoms[0].z, zmax = m_atoms[0].z;

    for (int i = 1; i < m_nAtoms; i++) {
    if (m_atoms[i].x < xmin) xmin = m_atoms[i].x;
    if (m_atoms[i].x > xmax) xmax = m_atoms[i].x;
    if (m_atoms[i].y < ymin) ymin = m_atoms[i].y;
    if (m_atoms[i].y > ymax) ymax = m_atoms[i].y;
    if (m_atoms[i].z < zmin) zmin = m_atoms[i].z;
    if (m_atoms[i].z > zmax) zmax = m_atoms[i].z;
    }

    // Add margin (EMmer uses 5A default)
    double margin = 5.0;
    xmin -= margin;
    xmax += margin;
    ymin -= margin;
    ymax += margin;
    zmin -= margin;
    zmax += margin;

    // Calculate grid dimensions
    m_nx = (int)ceil((xmax - xmin) / m_grid_spacing);
    m_ny = (int)ceil((ymax - ymin) / m_grid_spacing);
    m_nz = (int)ceil((zmax - zmin) / m_grid_spacing);

    // Ensure dimensions are at least 1
    if (m_nx < 1) m_nx = 1;
    if (m_ny < 1) m_ny = 1;
    if (m_nz < 1) m_nz = 1;

    m_origin[0] = xmin;
    m_origin[1] = ymin;
    m_origin[2] = zmin;
    }

    return 0;
}

double EmmerGenerator::calculateRefmacBlur() {
    // From GEMMI: blur = (8π²) / 1.1 * spacing² - b_min
    double spacing = m_d_min / (2 * m_config.rate);
    if (spacing <= 0) {
    spacing = m_grid_spacing;
    }

    // Find minimum B-factor
    double b_min = 1e10;
    for (int i = 0; i < m_nAtoms; i++) {
    if (m_atoms[i].b_factor < b_min) {
    b_min = m_atoms[i].b_factor;
    }
    }
    if (b_min < 0) b_min = 0;

    double u_to_b = 8.0 * M_PI * M_PI;
    double blur = (u_to_b / 1.1) * spacing * spacing - b_min;

    return (blur > 0) ? blur : 0.0;
}

void EmmerGenerator::precalculateForElement(EmmerElementData* elem, double B_eff) {
    const EmmerGaussianCoeff& coeff = COEFFICIENTS[elem->element_idx];

    const double _4pi = 4.0 * M_PI;

    for (int g = 0; g < EMMER_N_GAUSS; g++) {
    // t = 4π / (b + B)
    double t = _4pi / (coeff.b[g] + B_eff);

    // amplitude = a * t^(3/2)
    elem->precalc[g].amplitude = coeff.a[g] * pow(t, 1.5);

    // width = -t * π  (for exp(-width * r²))
    elem->precalc[g].width = -t * M_PI;

    // Estimate radius where density drops below cutoff
    // r = sqrt( log(amplitude/cutoff) / -width )
    double ratio = elem->precalc[g].amplitude / m_config.cutoff_level;
    if (ratio > 1.0) {
    elem->precalc[g].radius = sqrt(log(ratio) / -elem->precalc[g].width);
    }
    else {
    elem->precalc[g].radius = 0.0;
    }
    }
}

double EmmerGenerator::estimateRadius(const EmmerPrecalcGaussian* gaussians) {
    double max_radius = 0.0;
    for (int g = 0; g < EMMER_N_GAUSS; g++) {
    if (gaussians[g].radius > max_radius) {
    max_radius = gaussians[g].radius;
    }
    }
    return max_radius;
}

void EmmerGenerator::addAtomToGrid(int atom_idx, const EmmerElementData* elem, Ipp64f* grid) {
    const Atom& atom = m_atoms[atom_idx];

    // Convert to grid coordinates
    double gx = (atom.x - m_origin[0]) / m_grid_spacing;
    double gy = (atom.y - m_origin[1]) / m_grid_spacing;
    double gz = (atom.z - m_origin[2]) / m_grid_spacing;

    // Find integer coordinates
    int ix0 = (int)floor(gx);
    int iy0 = (int)floor(gy);
    int iz0 = (int)floor(gz);

    // Determine cutoff radius in voxels from the Gaussian with largest radius
    double max_radius = 0.0;
    for (int g = 0; g < EMMER_N_GAUSS; g++) {
        if (elem->precalc[g].radius > max_radius) {
            max_radius = elem->precalc[g].radius;
        }
    }
    int r_vox = (int)ceil(max_radius / m_grid_spacing) + 1;

    // Calculate bounds
    int imin = ix0 - r_vox;
    int imax = ix0 + r_vox;
    int jmin = iy0 - r_vox;
    int jmax = iy0 + r_vox;
    int kmin = iz0 - r_vox;
    int kmax = iz0 + r_vox;

    // Clamp to grid boundaries
    if (imin < 0) imin = 0;
    if (imax >= m_nx) imax = m_nx - 1;
    if (jmin < 0) jmin = 0;
    if (jmax >= m_ny) jmax = m_ny - 1;
    if (kmin < 0) kmin = 0;
    if (kmax >= m_nz) kmax = m_nz - 1;

    // Add Gaussian contribution to each voxel in the sphere
    for (int k = kmin; k <= kmax; k++) {
        double dz = (k - gz) * m_grid_spacing;
        double dz2 = dz * dz;

        for (int j = jmin; j <= jmax; j++) {
            double dy = (j - gy) * m_grid_spacing;
            double dy2 = dy * dy;

            int64_t base_idx = ((int64_t)k * m_ny + j) * m_nx;

            for (int i = imin; i <= imax; i++) {
                double dx = (i - gx) * m_grid_spacing;
                double r2 = dx * dx + dy2 + dz2;

                // Skip if outside the maximum radius
                if (r2 > max_radius * max_radius) continue;

                // Sum over 5 Gaussians
                double val = 0.0;
                for (int g = 0; g < EMMER_N_GAUSS; g++) {
                    val += elem->precalc[g].amplitude *
                        exp(elem->precalc[g].width * r2);
                }

                // Apply occupancy and add to grid
                int64_t idx = base_idx + i;
#pragma omp atomic
                grid[idx] += val * atom.occupancy;
            }
        }
    }
}

void EmmerGenerator::symmetrizeGrid() {
    // Simplified - in real implementation would use space group operations
    // For now, do nothing
}

void EmmerGenerator::applyOutputAlignment() {
    if (!m_config.align_output || m_nx <= 0 || m_ny <= 0 || m_nz <= 0) {
        return;
    }


    return;

    int64_t nvox = (int64_t)m_nx * m_ny * m_nz;

    // Create temporary grid for the transformed map
    Ipp64f* temp = ippsMalloc_64f(nvox);
    if (!temp) return;

    // Step 1: Flip along z-axis
    for (int k = 0; k < m_nz; k++) {
        int k_flip = m_nz - 1 - k;
        for (int j = 0; j < m_ny; j++) {
            for (int i = 0; i < m_nx; i++) {
                int64_t src_idx = ((int64_t)k * m_ny + j) * m_nx + i;
                int64_t dst_idx = ((int64_t)k_flip * m_ny + j) * m_nx + i;
                temp[dst_idx] = m_output[src_idx];
            }
        }
    }

    // Step 2: Rotate 90 degrees around y-axis (axes=(2,0) in EMmer)
    // New dimensions after rotation
    int new_nx = m_nz;
    int new_ny = m_ny;
    int new_nz = m_nx;

    int64_t new_nvox = (int64_t)new_nx * new_ny * new_nz;
    Ipp64f* rotated = ippsMalloc_64f(new_nvox);
    if (!rotated) {
        ippsFree(temp);
        return;
    }

    ippsZero_64f(rotated, new_nvox);

    // Correct rotation mapping: (x, y, z) -> (-z, y, x)
    for (int k = 0; k < m_nz; k++) {
        for (int j = 0; j < m_ny; j++) {
            for (int i = 0; i < m_nx; i++) {
                int64_t src_idx = ((int64_t)k * m_ny + j) * m_nx + i;
                if (temp[src_idx] == 0.0) continue;

                int new_i = m_nz - 1 - k;        // x' = -z (wrapped)
                int new_j = j;                     // y' = y
                int new_k = i;                      // z' = x

                if (new_i >= 0 && new_i < new_nx &&
                    new_j >= 0 && new_j < new_ny &&
                    new_k >= 0 && new_k < new_nz) {

                    int64_t dst_idx = ((int64_t)new_k * new_ny + new_j) * new_nx + new_i;
                    rotated[dst_idx] = temp[src_idx];
                }
            }
        }
    }

    // Calculate center of mass of atoms
    Ipp64f com_x = 0.0, com_y = 0.0, com_z = 0.0;
    int n_atoms = 0;

    for (int e = 0; e < m_nElements; e++) {
        for (int a = 0; a < m_element_data[e].atom_count; a++) {
            int atom_idx = m_element_data[e].atom_indices[a];
            com_x += m_atoms[atom_idx].x;
            com_y += m_atoms[atom_idx].y;
            com_z += m_atoms[atom_idx].z;
            n_atoms++;
        }
    }

    if (n_atoms > 0) {
        com_x /= n_atoms;
        com_y /= n_atoms;
        com_z /= n_atoms;
    }

    printf("\nOriginal center of mass: (%.3f, %.3f, %.3f)\n", com_x, com_y, com_z);

    // Transform center of mass: (x, y, z) -> (-z, y, x)
    Ipp64f new_com_x = -com_z;
    Ipp64f new_com_y = com_y;
    Ipp64f new_com_z = com_x;

    printf("New center of mass should be: (%.3f, %.3f, %.3f)\n",
        new_com_x, new_com_y, new_com_z);

    // Free old output and replace with rotated map
    ippsFree(m_output);
    m_output = rotated;

    // Update dimensions
    int old_nx = m_nx;
    int old_ny = m_ny;
    int old_nz = m_nz;

    m_nx = new_nx;
    m_ny = new_ny;
    m_nz = new_nz;

    // Calculate grid center in voxel coordinates
    Ipp64f grid_center_x = (m_nx - 1) * m_grid_spacing * 0.5;
    Ipp64f grid_center_y = (m_ny - 1) * m_grid_spacing * 0.5;
    Ipp64f grid_center_z = (m_nz - 1) * m_grid_spacing * 0.5;

    // Set origin so that the transformed center of mass is at the grid center
    m_origin[0] = new_com_x - grid_center_x;
    m_origin[1] = new_com_y - grid_center_y;
    m_origin[2] = new_com_z - grid_center_z;

    // Shift to ensure all coordinates are non-negative if needed
    Ipp64f min_allowed = 0.0;
    Ipp64f shift_x = (m_origin[0] < min_allowed) ? -m_origin[0] : 0.0;
    Ipp64f shift_y = (m_origin[1] < min_allowed) ? -m_origin[1] : 0.0;
    Ipp64f shift_z = (m_origin[2] < min_allowed) ? -m_origin[2] : 0.0;

    m_origin[0] += shift_x;
    m_origin[1] += shift_y;
    m_origin[2] += shift_z;

    printf("Final origin: (%.3f, %.3f, %.3f)\n",
        m_origin[0], m_origin[1], m_origin[2]);
    printf("Final center of mass position: (%.3f, %.3f, %.3f)\n",
        m_origin[0] + grid_center_x,
        m_origin[1] + grid_center_y,
        m_origin[2] + grid_center_z);
    printf("Final box bounds X: [%.3f, %.3f]\n",
        m_origin[0], m_origin[0] + (m_nx - 1) * m_grid_spacing);
    printf("Final box bounds Y: [%.3f, %.3f]\n",
        m_origin[1], m_origin[1] + (m_ny - 1) * m_grid_spacing);
    printf("Final box bounds Z: [%.3f, %.3f]\n",
        m_origin[2], m_origin[2] + (m_nz - 1) * m_grid_spacing);

    ippsFree(temp);
}

//=============================================================================
// Public methods
//=============================================================================

int EmmerGenerator::init(const EmmerGeneratorConfig* config,
    const Atom* atoms,
    int n_atoms,
    Ipp64f resolution,
    Ipp64f grid_spacing,
    const Ipp64f* cell_dimensions) {
    if (!config || !atoms || n_atoms <= 0 || resolution <= 0.0) {
    return -1;
    }

    m_config = *config;
    m_atoms = atoms;
    m_nAtoms = n_atoms;
    m_d_min = resolution;
    m_grid_spacing = (grid_spacing > 0) ? grid_spacing : resolution / 3.0;

    // Find unique elements
    int n_unique = findUniqueElements();
    if (n_unique <= 0) {
        return -2;
    }

    // Resolve grid parameters
    int ret = resolveParameters(cell_dimensions);
    if (ret != 0) {
        return -3;
    }

    // Calculate Refmac blur if requested
    if (m_config.set_refmac_blur && m_config.blur <= 0.0) {
        m_config.blur = calculateRefmacBlur();
    }

    // Allocate working grid for atom addition
    int64_t nvox = (int64_t)m_nx * m_ny * m_nz;
    m_work_grid = (Ipp64f*)ippsMalloc_64f(nvox);
    if (!m_work_grid) {
        return -4;
    }

    return 0;
}

int EmmerGenerator::run() {
    if (!m_atoms || m_nAtoms <= 0 || !m_work_grid) {
        return -1;
    }

    int64_t nvox = (int64_t)m_nx * m_ny * m_nz;

    // Allocate output grid
    m_output = (Ipp64f*)ippsMalloc_64f(nvox);
    if (!m_output) {
        return -2;
    }

    ippsZero_64f(m_output, nvox);

    printf("\n=== EMmer-style Map Generation ===\n");
    printf("Resolution: %.2f A\n", m_d_min);
    printf("Grid spacing: %.2f A\n", m_grid_spacing);
    printf("Grid dimensions: %d x %d x %d (%lld voxels)\n", m_nx, m_ny, m_nz, nvox);
    printf("Refmac blur: %.2f A^2\n", m_config.blur);
    printf("Atoms: %d\n", m_nAtoms);
    printf("Unique elements: %d\n", m_nElements);

    // Precalculate Gaussians for each element type
    for (int e = 0; e < m_nElements; e++) {
        EmmerElementData* elem = &m_element_data[e];

        // Calculate average B-factor for this element
        double avg_B = 0.0;
        for (int a = 0; a < elem->atom_count; a++) {
            int atom_idx = elem->atom_indices[a];
            avg_B += m_atoms[atom_idx].b_factor;
        }
        if (elem->atom_count > 0) {
            avg_B /= elem->atom_count;
        }
        elem->avg_b_factor = avg_B;

        precalculateForElement(elem, avg_B + m_config.blur);

        printf("  Element %s (Z=%d): %d atoms, avg B=%.2f\n",
        COEFFICIENTS[elem->element_idx].element,
        COEFFICIENTS[elem->element_idx].atomic_number,
        elem->atom_count, avg_B);
    }

    // Add atoms to grid
    printf("\nAdding atoms to grid...\n");

    for (int e = 0; e < m_nElements; e++) {
        EmmerElementData* elem = &m_element_data[e];

#pragma omp parallel for
        for (int a = 0; a < elem->atom_count; a++) {
            int atom_idx = elem->atom_indices[a];
            addAtomToGrid(atom_idx, elem, m_output);
        }
    }

    // Apply symmetry if requested
    if (m_config.symmetry_expansion) {
        symmetrizeGrid();
    }

    // Apply output alignment if requested
    if (m_config.align_output) {
        applyOutputAlignment();
    }

    // Calculate statistics
    Ipp64f min_val, max_val, mean_val, std_val;
    ippsMinMax_64f(m_output, nvox, &min_val, &max_val);
    ippsMean_64f(m_output, nvox, &mean_val);
    ippsStdDev_64f(m_output, nvox, &std_val);

    printf("\nMap statistics:\n");
    printf("  Min: %.6f\n", min_val);
    printf("  Max: %.6f\n", max_val);
    printf("  Mean: %.6f\n", mean_val);
    printf("  Std Dev: %.6f\n", std_val);
    printf("================================\n");

    return 0;
}

void EmmerGenerator::getOutputOrigin(Ipp64f origin[3]) const {
    origin[0] = m_origin[0];
    origin[1] = m_origin[1];
    origin[2] = m_origin[2];
}

Ipp64f* EmmerGenerator::releaseOutput() {
    Ipp64f* map = m_output;
    m_output = nullptr;
    return map;
}