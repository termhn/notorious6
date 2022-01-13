// TL;DR: Compress brightness, preserve chroma by desaturation.

#include "inc/prelude.glsl"
#include "inc/ictcp.hlsl"
#include "inc/luv.hlsl"
#include "inc/lms.hlsl"
#include "inc/oklab.hlsl"
#include "inc/lab.hlsl"
#include "inc/h-k.hlsl"
#include "inc/ycbcr.hlsl"

// Helmholtz-Kohlrausch adjustment methods
#define HK_ADJUSTMENT_METHOD_NAYATANI 0
#define HK_ADJUSTMENT_METHOD_NONE 1

// Brightness compression curves:
#define BRIGHTNESS_COMPRESSION_CURVE_REINHARD 0
#define BRIGHTNESS_COMPRESSION_CURVE_SIRAGUSANO_SMITH 1    // :P

// ----------------------------------------------------------------
// Configurable stuff:

#define BRIGHTNESS_COMPRESSION_CURVE BRIGHTNESS_COMPRESSION_CURVE_SIRAGUSANO_SMITH

// Choose the method for performing the H-K adjustment
#define HK_ADJUSTMENT_METHOD HK_ADJUSTMENT_METHOD_NAYATANI

// Adapting luminance (L_a) used for the H-K adjustment. 20 cd/m2 was used in Sanders and Wyszecki (1964)
#define HK_ADAPTING_LUMINANCE 20

// Match target compressed brightness while attenuating chroma.
// Important in the low end, as well as at the high end of blue and red.
#define USE_BRIGHTNESS_LINEAR_CHROMA_ATTENUATION 1

// The amount of the outer part of the output gamut to use for gamut chroma compression.
// For example, 0.25 will use the outer 25% of the gamut to compress out-of gamut + in-gamut colors into.
// Thus, you'll lose some chroma for colors which are inside the original part of this section of the output gamma
// so that out of gamma colors may be compressed into it while retaining relational stability.
#define GAMUT_COMPRESS_CUSP 0.2 // compression at max chroma for that hue
#define GAMUT_COMPRESS_SHADOW 0.05 // compression at 1.0 lightness
#define GAMUT_COMPRESS_HIGHLIGHT 0.5 // compression at 0 lightness

// A higher power will result in biasing higher compression towards the outside of the output gamut,
// therefore higher compression ratios for colors further outside
// the output gamut.
#define GAMUT_COMPRESS_POWER 2.0

// Should stay between > 0. Shifts the color compression curve similar to the power, but affects the more inner
// parts of the gamut being compressed.
#define GAMUT_COMPRESS_BIAS 0.1

// Controls for manual desaturation of lighter than "white" stimulus (greens, yellows);
// see comments in the code for more details.
#define GLOBAL_CHROMA_ATTENUATION_START 0.0
#define GLOBAL_CHROMA_ATTENUATION_EXPONENT 4.0
// ----------------------------------------------------------------


// Map (hk-adjusted) scene-linear luminance through a curve yielding values in 0..1, yielding display-linear luminance
float compress_luminance(float v) {
	#if BRIGHTNESS_COMPRESSION_CURVE == BRIGHTNESS_COMPRESSION_CURVE_REINHARD
		// Reinhard
		return v / (v + 1.0);
	#elif BRIGHTNESS_COMPRESSION_CURVE == BRIGHTNESS_COMPRESSION_CURVE_SIRAGUSANO_SMITH
		// From Jed Smith: https://github.com/jedypod/open-display-transform/wiki/tech_tonescale,
        // based on stuff from Daniele Siragusano: https://community.acescentral.com/t/output-transform-tone-scale/3498/14
        // Reinhard with flare compensation.
        const float scene_exposure = 1.0;
        const float toe_power = 1.2;
        const float display_exposure = 1.0205;
		return saturate(display_exposure * pow(v / (v + scene_exposure), toe_power));
    #endif
}

float linear_srgb_to_luminance(float3 col) {
    return rgb_to_ycbcr(col).x;
}

// Scene-linear luminance adjusted by the Helmholtz-Kohlrausch effect
float srgb_to_hk_adjusted_luminance(float3 input) {
#if HK_ADJUSTMENT_METHOD == HK_ADJUSTMENT_METHOD_NAYATANI
    const float luminance = linear_srgb_to_luminance(input);
    const float2 uv = cie_XYZ_to_Luv_uv(RGBToXYZ(input));
    const float luv_brightness = hsluv_yToL(luminance);
    const float mult = nayatani_hk_lightness_adjustment_multiplier(uv, HK_ADAPTING_LUMINANCE);
    return hsluv_lToY(luv_brightness * mult);
#elif HK_ADJUSTMENT_METHOD == HK_ADJUSTMENT_METHOD_NONE
    return linear_srgb_to_luminance(input);
#endif
}

// Based on https://bottosson.github.io/posts/gamutclipping/
float3 gamut_compress_to_achromatic_luminance(float3 rgb, float achromatic_luminance)
{
	float3 lab = linear_srgb_to_oklab(rgb);

    // The achromatic stimulus we'll interpolate towards to fix out-of-gamut stimulus.
    float L0 = saturate(linear_srgb_luminance_to_oklab_L(achromatic_luminance));

    // Values lighter than "white" are already within the gamut, so our brightness compression is "done".
    // Perceptually they look wrong though, as they don't follow the desaturation that other stimulus does.
    // We fix that manually here by biasing the interpolation towards "white" at the end of the brightness range.
    // This "fixes" the yellows and greens.
    const float supplemental_chroma_attenuation = pow(
        saturate(
            (achromatic_luminance - 1.0 * GLOBAL_CHROMA_ATTENUATION_START)
            / (1.0 - GLOBAL_CHROMA_ATTENUATION_START)
        ), GLOBAL_CHROMA_ATTENUATION_EXPONENT
    );

    lab = lerp(lab, float3(L0, 0.0, 0.0), supplemental_chroma_attenuation);

	float L = lab.x;
	float C = length(lab.yz);
    // if already achromatic, just clip and return
    if (C <= 0.00001) {
        // return float3(1.0, 0.0, 1.0);
        return saturate(achromatic_luminance).xxx;
    }
	float a_norm = lab.y / C;
	float b_norm = lab.z / C;

	// Find the cusp of the gamut triangle
	float2 cusp = find_cusp(a_norm, b_norm);
    float L_cusp = cusp.x;
    float C_cusp = cusp.y;

	float t_gamut_boundary = find_gamut_intersection(a_norm, b_norm, cusp, L, C, L0);

    // gamut boundary at 0, clip and return achromatic
    if (t_gamut_boundary <= 0.0) {
        return saturate(achromatic_luminance).xxx;
    }

    float blend_to_top = (L0 - L_cusp) / (1.0 - L_cusp);
    float blend_to_bottom = 1.0 - (L0 / L_cusp);

    float blend = L0 > L_cusp ? blend_to_top : blend_to_bottom;
    float compress_at_pole = L0 > L_cusp ? GAMUT_COMPRESS_HIGHLIGHT : GAMUT_COMPRESS_SHADOW;
    float compress = lerp(GAMUT_COMPRESS_CUSP, compress_at_pole, blend);

    // return oklab_to_linear_srgb(float3(compress, 0.0, 0.0));
    float t_compression_boundary = t_gamut_boundary * (1.0 - compress);

    float vector_len = length(float2(L, C));

    float delta_len = max(0.0, (1.0 - t_compression_boundary) * vector_len);
    
    // adjusted from https://www.desmos.com/calculator/54aytu7hek
    float compression = max(0.0, delta_len / pow(GAMUT_COMPRESS_BIAS + pow(delta_len, GAMUT_COMPRESS_POWER), 1/GAMUT_COMPRESS_POWER));

    float t = 1.0;
    if (t_compression_boundary < 1.0) {
        t = lerp(t_compression_boundary, t_gamut_boundary, compression);
    }

    // float vector_len = length(float2(L, C));

    // remap so that 1.0 = gamut boundary
    // float x = 1.0 / (vector_len * t_gamut_boundary);

    // float compressed = 0.2;
    // float protected = 1.0 - compressed;
    // float p = 2.5;
    // float x_t_s = (x - protected) / compressed;

    // float t = x;
    // if (x > protected) {
    //     t = protected + compressed * x_t_s / pow(pow(x_t_s, p) + 1, 1 / p);
    // }

    // t = t * t_gamut_boundary * vector_len;

	float L_comp = lerp(L0, L, t);
    float C_comp = C * t;
	return oklab_to_linear_srgb(float3(L_comp, C_comp * a_norm, C_comp * b_norm));
}


float3 compress_stimulus(float3 input) {
    // Find the input luminance adjusted by the Helmholtz-Kohlrausch effect.
    const float input_luminance_hk_adjusted = srgb_to_hk_adjusted_luminance(input);

    // The highest displayable intensity stimulus with the same chromaticity as the input,
    // and its associated HK-adjusted luminance.
    const float3 max_displayable_rgb = input / max(input.r, max(input.g, input.b)).xxx;
    float max_displayable_rgb_lum_hk_adjusted = srgb_to_hk_adjusted_luminance(max_displayable_rgb);

    // Compress the luminance. We will then adjust the chromatic input stimulus to match this brightness.
    // Note that this is not the non-linear "L*", but a 0..1 value as a multilpier
    // over the maximum achr(hk-adjusted) luminanceinance.
    const float compressed_achromatic_luminance = compress_luminance(input_luminance_hk_adjusted);

    // Scale the chromatic stimulus so that its luminance matches `compressed_achromatic_luminance`.
    // TODO: Overly simplistic, and does not accurately map the brightness.
    //
    // This will create (mostly) matching brightness, but potentially out of gamut components.
    float3 compressed_rgb = (max_displayable_rgb / max_displayable_rgb_lum_hk_adjusted) * compressed_achromatic_luminance;


    // We now want to map the out-of-gamut stimulus back to what our device can display.
    // Since both the `compressed_rgb` and `clamped_compressed_achromatic_luminance` are of the same-ish
    // brightness, and `clamped_compressed_achromatic_luminance.xxx` is guaranteed to be inside the gamut,
    // we can trace a path from `compressed_rgb` towards `clamped_compressed_achromatic_luminance.xxx`,
    // and stop once we have intersected the target gamut.

    // This has the effect of removing chromatic content from the compressed stimulus,
    // and replacing that with achromatic content. If we do that naively, we run into
    // a perceptual hue shift due to the Abney effect.
    //
    // To counter, we first transform both vertices of the path we want to trace
    // into a perceptual space which preserves sensation of hue, then we trace
    // a straight line _inside that space_ until we intersect the gamut.

    compressed_rgb = gamut_compress_to_achromatic_luminance(compressed_rgb, compressed_achromatic_luminance);

    //return srgb_to_hk_adjusted_luminance(compressed_rgb).xxx;
    return compressed_rgb;
}
