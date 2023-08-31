extern crate nalgebra as na;
// use ndarray::{Array1};
// use crate::fftconv_externo::fftconvolve;

use na::{DVector,Matrix2};

use fft_convolver::FFTConvolver;

// pub fn fftconvolve_envuelto(
//     in1: Vec<f32>,
//     in2: Vec<f32>,
// ) -> Vec<f32> 
// {
//     let in1 = Array1::<f32>::from_vec(in1);
//     let in2 = Array1::<f32>::from_vec(in2);
//     fftconvolve(&in1, &in2).unwrap().to_vec()
// }

pub fn fftconvolve_envuelto2(
    impulse_response: Vec<f32>,
    input: Vec<f32>,
) -> Vec<f32> 
{
    let mut convolver = FFTConvolver::default();
    convolver.init( 16 , &impulse_response);
    let mut output = vec![0_f32 ; input.len()];

    convolver.process(&input, &mut output);
    output
}

//
// Convierte la primera fila de una matriz de covarianza toeplitz (U) en 
// la función de impulso que genera el correspondiente proceso gaussiano,
// que casualmente es la última fila de la L cuando hacemos cholesky para
// descomponer U en LL^T con L lower-triangular.
//
pub fn toeplitz_covariance_to_impulse(
    cov: Vec<f32>,
    ) -> Vec<f32>
{

    let m = cov.len();
    let mut uZu: DVector<f32> = DVector::repeat(m+1,0.0);

    {
        let mut u = uZu.rows_mut(1,m);
        u /= u[0].sqrt();
    }

    let mut v: DVector<f32> = uZu.rows(0,m).clone_owned();

    let mut h: DVector<f32> = DVector::repeat(m,0.0);
    h[m-1] = 0.;

    for k in 0..m-1
    {

        let s: f32 = v[k+1]/uZu[k+1];
        let c: f32 = (1.0-f32::powi(s,2)).sqrt();

        let mut new_u = (uZu.rows(0,m) - s * v.clone()) / c;

        {
        let mut u = uZu.rows_mut(1,m);
        u = new_u.rows_mut(0,m);
        v = (- s * u) + c * v;
        }
        
        h[m-1-k] = uZu[m]
    }

    h.data.into()
}

// create tests
#[cfg(test)]
mod tests {
    use super::*;
    // use ndarray::{Array1, Array2};

    // #[test]
    // fn test_fftconvolve_1d() {
    //     let in1: Vec<f32> = vec![1.0, 2.0, 3.0];
    //     let in2: Vec<f32> = vec![6.0, 5.0, 4.0];
    //     let out = fftconvolve_envuelto(in1, in2);

    //     let expected: Vec<f32> = vec![6., 17., 32., 23., 12.];

    //     assert_eq_float_vec!(out, expected, 0.000001);

    // }

    #[test]
    fn test_fftconvolve2_1d() {
        let in1: Vec<f32> = vec![1.0, 2.0, 3.0];
        let in2: Vec<f32> = vec![6.0, 5.0, 4.0];
        let out = fftconvolve_envuelto2(in1, in2);

        // El modo de esta convolución no es "valid", no importa.
        let expected: Vec<f32> = vec![6., 17., 32.];

        assert_eq_float_vec!(out, expected, 0.00001);

    }
    #[test]
    fn test_toeplitz_covariance_to_impulse() {
        let input: Vec<f32> = vec![1.0,0.1,0.2,0.3,0.4];
        let out = toeplitz_covariance_to_impulse(input);

        let expected: Vec<f32> = vec![ 0.        ,  0.86964539, -0.07181627,  0.10137304,  0.26130983];

        assert_eq_float_vec!(out, expected, 0.00001);
    }
}
