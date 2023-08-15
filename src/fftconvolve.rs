use ndarray::{Array1};
use crate::fftconv_externo::fftconvolve;

pub fn fftconvolve_envuelto(
    in1: Vec<f32>,
    in2: Vec<f32>,
) -> Vec<f32> 
{
    let in1 = Array1::<f32>::from_vec(in1);
    let in2 = Array1::<f32>::from_vec(in2);
    fftconvolve(&in1, &in2).unwrap().to_vec()
}

// create tests
#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::{Array1, Array2};

    #[test]
    fn test_fftconvolve_1d() {
        let in1: Vec<f32> = vec![1.0, 2.0, 3.0];
        let in2: Vec<f32> = vec![6.0, 5.0, 4.0];
        let out = fftconvolve_envuelto(in1, in2);

        let expected: Vec<f32> = vec![6., 17., 32., 23., 12.];

        assert_eq_float_vec!(out, expected, 0.000001);

    }
}
