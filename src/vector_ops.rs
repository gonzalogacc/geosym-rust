extern crate nalgebra as na;

use na::{DVector};

use fft_convolver::FFTConvolver;

pub fn fftconvolve_envuelto(
    impulse_response: Vec<f32>,
    input: Vec<f32>,
) -> Vec<f32> 
{
    let mut convolver = FFTConvolver::default();

    // Acá podría haber error, manejarlo!
    let _ = convolver.init( 16 , &impulse_response);
    let mut output = vec![0_f32 ; input.len()];

    // Acá podría haber error, manejarlo!
    let _ =convolver.process(&input, &mut output);

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

    // Quiero un 0 adelante de u, porque a veces quiero mirar ese cero, que en el paper es Z u y
    // aveces no, lo que en el paper es u.
    let mut zu: DVector<f32> = DVector::repeat(m+1,0.0);

    {
        let mut u = zu.rows_mut(1,m);
        u += DVector::from_vec(cov);
        u /= u[0].sqrt();
    }

    // 
    let mut v: DVector<f32> = zu.rows(1,m).into_owned();
    v[0] = 0.0;

    let mut h: DVector<f32> = DVector::repeat(m,0.0);
    h[m-1] = 0.;

    for k in 0..m-1
    {

        // Textual del paper, aunque u_k es Zu_k+1
        let s: f32 = v[k+1]/zu[k+1];
        let c: f32 = (1.0-f32::powi(s,2)).sqrt();

        // Esta es la cuenta del paper Zu - seno...
        let new_u = (zu.rows(0,m) - s * v.clone()) / c;

        {
        // Ridicula forma de actualizar los valores de u
        // porque no entiendo como asignar a un slice
        // de nalgebra.
        let mut u = zu.rows_mut(1,m);
        u *= 0.0;
        u += new_u.rows(0,m);
        v = (- s * u) + c * v;
        }

        h[m-1-k] = zu[m]
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

    fn toeplitz_matrix( v: Vec<f32> ) -> na::DMatrix<f32>
    {
        let m = v.len();
        //let largo=np.concatenate([v[::-1],v[1:]])
        let largo: Vec<f32> = v.clone().into_iter().rev().chain(
                              v.clone()[1..m].into_iter().copied()).collect();

        let mx: Vec<f32> = (0..m).map(|i| largo.clone()[m-1-i..2*m-i-1]
                                            .into_iter().copied().collect::<Vec<f32>>())
                                   .flatten().collect();
        na::DMatrix::from_vec(m,m,mx)
    }

    fn cholesky_toeplitz_into_impulse( v: Vec<f32> ) -> Vec<f32>
    {
        let m = v.len();
        let t = toeplitz_matrix(v);
        //let t_data : Vec<f32> = t.cholesky().expect("La Matriz no admite descomposición LU")
        //                         .l().rows(m-1,1).into_owned().data.as_slice().iter()
        //                         .rev().copied().collect();
        let t_data : Vec<f32> = t.cholesky().expect("La Matriz no admite descomposición LU")
                                   .l().rows(m-1,1).into_owned().data.as_slice().iter()
                                   .chain(&[0.0])
                                   .rev().take(m).copied().collect();
        t_data
    }

    #[test]
    fn test_fftconvolve2_1d() {
        let in1: Vec<f32> = vec![1.0, 2.0, 3.0];
        let in2: Vec<f32> = vec![6.0, 5.0, 4.0];
        let out = fftconvolve_envuelto(in1, in2);

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

    #[test]
    fn test_toeplitz_covariance_to_impulse_cholesky() {
        //
        // Cualquier vector que empiece con un numero y tenga numeros menores
        // después, ordenados de menor a mayor, forma una matriz valida.
        
        let input: Vec<f32> = vec![1.0,0.1,0.2,0.3,0.4];
        let out = toeplitz_covariance_to_impulse(input.clone());
        let expected: Vec<f32> = cholesky_toeplitz_into_impulse(input);

        assert_eq_float_vec!(out, expected, 0.00001);

        let input: Vec<f32> = vec![3.8,0.1,0.2,0.35,0.4,0.42,0.48];
        let out = toeplitz_covariance_to_impulse(input.clone());
        let expected: Vec<f32> = cholesky_toeplitz_into_impulse(input);

        assert_eq_float_vec!(out, expected, 0.00001);

    }
}
