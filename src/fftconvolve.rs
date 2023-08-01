use easyfft::prelude::*;
use easyfft::{Complex, FftNum};
use ndarray::{OwnedRepr};
use ndarray::{Array1,Array, ArrayBase, Axis, Data, Slice, Ix1};
use num::FromPrimitive;
use num_traits::Zero;
use std::error::Error;

pub fn fftconvolve_vec(
    in1: Vec<f32>,
    in2: Vec<f32>,
) -> Vec<f32> 
{
    let in1 = Array1::<f32>::from_vec(in1);
    let in2 = Array1::<f32>::from_vec(in2);
    fftconvolve(&in1, &in2).unwrap().to_vec()
}

// Determine the size of the output of convolution
// Only full mode implemented
// pub enum Mode {
//     Full,
//     Same,
//     Valid,
// }

// Pad the edges of an array with zeros.
//
// `pad_width` specifies the length of the padding at the beginning
// and end of each axis.
//
// Returns a Result. Errors if arr.ndim() != pad_width.len().
fn pad_with_zeros<A, S>(
    arr: &ArrayBase<S, Ix1>,
    pad_width: [usize; 2],
) -> Result<Array<A, Ix1>, Box<dyn Error>>
where
    A: FftNum,
    S: Data<Elem = A>,
{
    if arr.ndim() > 1 {
        return Err("arr.ndim() > 1".into());
    }

    // Compute shape of final padded array.
    let mut padded_shape = arr.raw_dim();
    {
        let ax = 0;
        let ax_len = arr.shape()[0];
        let [pad_lo, pad_hi] = pad_width;

        padded_shape[ax] = ax_len + pad_lo + pad_hi;
    } 

    let mut padded = Array::zeros(padded_shape);
    let padded_dim = padded.raw_dim();
    {
        // Select portion of padded array that needs to be copied from the
        // original array.
        let mut orig_portion = padded.view_mut();
        {
            let ax = 0;
            let [pad_lo, pad_hi] = pad_width;

            orig_portion.slice_axis_inplace(
                Axis(ax),
                Slice::from(pad_lo as isize..padded_dim[ax] as isize - (pad_hi as isize)),
            );
        }
        // Copy the data from the original array.
        orig_portion.assign(arr);
    }
    Ok(padded)
}

/// Generates a Vec<[usize; 2]> specifying how much to pad each axis.
fn generate_pad_vector<A, S>(arr: &ArrayBase<S, Ix1>, shape: usize) -> [usize; 2]
where
    A: FftNum,
    S: Data<Elem = A>,
{
    let arr_size = arr.shape()[0];
    let new_size = shape;
    [0, new_size - arr_size]
}

/// Convolve two N-dimensional arrays using FFT.
fn fftconvolve<A, S>(
    in1: &ArrayBase<S, Ix1>,
    in2: &ArrayBase<S, Ix1>,
) -> Result<ArrayBase<OwnedRepr<A>, Ix1>, Box<dyn Error>> 
where
    A: FftNum + FromPrimitive,
    S: Data<Elem = A>,
{
    // check that arrays have the same number of dimensions
    if in1.ndim() != in2.ndim() {
        return Err("Input arrays must have the same number of dimensions.".into());
    }

    // // Pad the arrays to the next power of 2.
    let mut shape = in1.shape()[0];

    shape = shape + in2.shape()[0] - 1;

    let in1 = pad_with_zeros(in1, generate_pad_vector(in1, shape))?;
    let in2 = pad_with_zeros(in2, generate_pad_vector(in2, shape))?;

    // multiple values in shape together to get total size
    let total_size = shape;

    // Crea un array con complejos para cada input
    let mut in1 = in1.mapv(|x| Complex::new(x, Zero::zero()));
    let mut in2 = in2.mapv(|x| Complex::new(x, Zero::zero()));

    // fft_mut es parte del trait FftNum, que está implementado
    // sobre slices de complejos en easyfft.
    // Hace la dfft in-place.
    in1.as_slice_mut().unwrap().fft_mut();
    in2.as_slice_mut().unwrap().fft_mut();

    // Multiply the FFT.
    let mut out = in1 * in2;

    // Lo mismo, la ifft está implementada para slices de complejos en el
    // trait FftNum.
    out.as_slice_mut().unwrap().ifft_mut();

    // Return the real part of the result. Note normalise by 1/total_size
    let total_size = A::from_usize(total_size).unwrap();

    let out = out.mapv(|x| x.re / total_size);
    Ok(out)
}

// create tests
#[cfg(test)]
mod tests {
    use super::*;
    use ndarray_linalg::assert_aclose;
    use ndarray::{Array1, Array2};

    #[test]
    fn reverse_array() {
        let mut standard = Array2::from_shape_vec((3, 3), vec![1, 2, 3, 4, 5, 6, 7, 8, 9]).unwrap();
        // let mut standard = standard.into_iter().collect::<Vec<_>>();
        // standard.reverse();
        standard.slice_each_axis_inplace(|_| Slice::new(0, None, -1));
        // reverse axes
        // let reversed = Array2::from_shape_vec((3, 3), standard).unwrap();
        let expected = Array2::from_shape_vec((3, 3), vec![9, 8, 7, 6, 5, 4, 3, 2, 1]).unwrap();
        assert_eq!(standard, expected);
    }

    #[test]
    fn test_pad() {
        let to_pad = Array1::<f32>::ones((4,));
        let padded = pad_with_zeros(&to_pad, [3, 3]).unwrap();
        let expected = Array1::<f32>::from_shape_vec(
            (10,),
            vec![
                0., 0., 0., 1., 1., 1., 1., 0., 0., 0.
            ],
        )
        .unwrap();
        // assert that the padded array is equal to the expected array
        assert_eq!(padded, expected);
    }

    #[test]
    fn test_fftconvolve_1d() {
        let in1 = Array1::range(1.0, 4.0, 1.0);
        let in2 = Array1::range(6.0, 3.0, -1.0);
        let out = fftconvolve(&in1, &in2).unwrap();
        let expected = Array1::<f64>::from_vec(vec![6., 17., 32., 23., 12.]);
        out.iter()
            .zip(expected.iter())
            .for_each(|(a, b)| assert_aclose!(*a, *b, 1e-6));
    }

    #[test]
    fn test_fftconvolve_1d_32() {
        let in1 = Array1::<f32>::range(1.0, 4.0, 1.0);
        let in2 = Array1::<f32>::range(6.0, 3.0, -1.0);
        let out = fftconvolve(&in1, &in2).unwrap();
        let expected = Array1::<f32>::from_vec(vec![6., 17., 32., 23., 12.]);
        out.iter()
            .zip(expected.iter())
            .for_each(|(a, b)| assert_aclose!(*a, *b, 1e-6));
    }
}
