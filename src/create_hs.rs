use std::io;
use std::error::Error;

pub struct SigmaH {
    pub sigma: f32,
    pub h: Vec<f32>,
}

impl SigmaH {
    fn new(sigma: f32, h: Vec<f32>) -> SigmaH {
        SigmaH {
            sigma,
            h,
        }
    }
}

pub fn white_noise(m: usize, sigma: f32) -> Result<SigmaH, Box<dyn Error>> {
    // Create impulse function for White noise.
    let mut h: Vec<f32> = vec![0.0; m];
    h[0] = sigma;
    Ok(SigmaH::new(sigma, h))
}

fn recursion_power_flicker_rw(m: usize, d: f32) -> Vec<f32>{
    // Recursion to create impulse function for Powerlay, Flicker or RW noise
    // Flicker is Powerlaw with spectral density 0.5
    // RandomWalk is Powerlaw with spectral density 1.0
    let mut h: Vec<f32> = vec![0.0; m];
    
    h[0] = 1.0;
    let mut h0: f32 = 1.0;
    for i in 1..m {
        h[i] = (d+i as f32-1.0)/i as f32 * h0;
        h0 = (d+i as f32-1.0)/i as f32 * h0;
    }
    h
}

macro_rules! assert_eq_float {
    ($v1: expr, $v2: expr, $d: expr ) => {
        if ($v1 - $v2).abs() > $d {
            panic!("Values {:?} != {:?} differ more than {:?}", $v1, $v2, $d);
        }
    }
}

macro_rules! assert_eq_float_vec {
    ($v1: expr, $v2: expr, $d: expr ) => {
        if $v1.len() != $v2.len() {
            panic!("Vectors have different length {:?} : {:?}", $v1.len(), $v2.len());
        }
        for (i, value) in $v1.iter().enumerate(){
            println!("{:?}, {:?}", i, value);
            if (value - $v2[i]).abs() > $d {
                panic!("Values {:?} != {:?} differ more than {:?}", value, $v2[i], $d);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::create_hs::{white_noise, recursion_power_flicker_rw};

    #[test]
    fn test_white_ok() {
        let sigma = 1.0;
        let m: usize = 10;
        let sigma_h = white_noise(m, sigma).unwrap();
        assert_eq!(sigma_h.sigma, sigma);
        assert_eq!(sigma_h.h, vec![sigma, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    }
    
    #[test]
    fn test_recursion_power_flicker_rw_ok() {
        // Write other tests
        let d: f32 = 0.5;
        let m: usize = 10;
        let response = recursion_power_flicker_rw(m, d);
        // Taken from the python function (run the function in a shell)
        let expected: Vec<f32> = vec![1., 0.5, 0.375 , 0.3125, 0.2734375 ,0.24609375, 0.22558594, 0.20947266, 0.19638062, 0.18547058];
        assert_eq!(response, expected);
    }

    #[test]
    fn test_recursion_power_flicker_rw_100() {
        // Write other tests
        let d: f32 = 0.1;
        let m: usize = 100;
        let response = recursion_power_flicker_rw(m, d);
        
        // Taken from the python function (run the function in a shell)
        let expected: Vec<f32> = vec![1., 0.1 ,0.055 ,0.0385 ,0.0298375 ,
       0.02446675, 0.02079674, 0.01812287, 0.01608405, 0.01447564,
       0.01317284, 0.01209506, 0.01118793, 0.01041338, 0.00974395,
       0.00915931, 0.0086441 , 0.00818647, 0.00777715, 0.00740876,
       0.00707536, 0.00677213, 0.00649509, 0.00624094, 0.0060069 ,
       0.00579065, 0.00559021, 0.00540387, 0.00523017, 0.00506785,
       0.00491582, 0.0047731 , 0.00463886, 0.00451234, 0.0043929 ,
       0.00427994, 0.00417294, 0.00407144, 0.00397501, 0.00388328,
       0.0037959 , 0.00371258, 0.00363302, 0.00355698, 0.00348423,
       0.00341454, 0.00334774, 0.00328363, 0.00322206, 0.00316288,
       0.00310595, 0.00305114, 0.00299833, 0.00294742, 0.00289829,
       0.00285087, 0.00280505, 0.00276076, 0.00271792, 0.00267646,
       0.00263631, 0.00259742, 0.00255971, 0.00252314, 0.00248766,
       0.00245322, 0.00241976, 0.00238726, 0.00235566, 0.00232494,
       0.00229505, 0.00226595, 0.00223763, 0.00221004, 0.00218316,
       0.00215697, 0.00213142, 0.00210651, 0.0020822 , 0.00205848,
       0.00203532, 0.00201271, 0.00199062, 0.00196903, 0.00194794,
       0.00192731, 0.00190714, 0.00188741, 0.00186811, 0.00184922,
       0.00183073, 0.00181262, 0.00179489, 0.00177752, 0.0017605 ,
       0.00174382, 0.00172747, 0.00171145, 0.00169573, 0.00168031];
         
        assert_eq_float_vec!(response, expected, 0.000001);

    }
}
