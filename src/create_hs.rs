
use std::io;
use std::f32::consts::PI;
//use std::error::Error;
//use crate::utils::Units;
use serde::Deserialize;
use crate::vector_ops::toeplitz_covariance_to_impulse;

#[allow(dead_code)]
#[derive(Debug)]
pub enum ModelError {
    IOError(io::Error),
}

#[allow(dead_code)]
#[derive(Deserialize, Debug, Copy, Clone)]
pub enum Units {
    Mom,
    Msf,
}

#[derive(Debug)]
pub struct SigmaH {
    pub sigma: f32,
    pub h: Vec<f32>,
}

pub fn white(m: usize, sigma: f32) -> Result<SigmaH, ModelError> {
    // Create impulse function for White noise.
    let mut h: Vec<f32> = vec![0.0; m];
    // Tomamos la convención de que la función de impulso se escala antes de la convolución.
    // h[0] = sigma; 
    h[0] = 1.0;
    Ok(SigmaH {sigma, h})
}

pub fn powerlaw(m: usize, kappa: f32, sigma: f32, units: Units, dt: f32) -> Result<SigmaH, ModelError> {
    let d = -kappa/2.0;
    let gmsv = gauss_markov_scale_variance(sigma, d, units, dt).unwrap();
    let rpf = recursion_power_flicker_rw(m, d).unwrap();
    Ok(SigmaH {sigma: gmsv, h: rpf})
}

pub fn flicker(m: usize, sigma: f32, units: Units, dt: f32) -> Result<SigmaH, ModelError> {
    let d = 0.5;
    let gmsv = gauss_markov_scale_variance(sigma, d, units, dt).unwrap();
    let rpf = recursion_power_flicker_rw(m, d).unwrap();
    Ok(SigmaH {sigma: gmsv, h: rpf})
}

pub fn random_walk(m: usize, sigma: f32, units: Units, dt: f32) -> Result<SigmaH, ModelError> {
    let d = 1.0;
    let gmsv = gauss_markov_scale_variance(sigma, d, units, dt).unwrap();
    let rpf = recursion_power_flicker_rw(m, d).unwrap();
    Ok(SigmaH {sigma: gmsv, h: rpf})
}

pub fn generalized_gauss_markov(m: usize, kappa: f32, one_minus_phi: f32, sigma: f32, units: Units, dt: f32) -> Result<SigmaH, ModelError> {
    let d = -kappa/2.0;
    let gmsv = gauss_markov_scale_variance(sigma, d, units, dt).unwrap();
    let rpf = recursion_ggm(m, d, one_minus_phi).unwrap();
    Ok(SigmaH {sigma: gmsv, h: rpf})
}

pub fn varying_anual(m: usize, phi: f32, units: Units, dt: f32, sigma: f32) -> Result<SigmaH, ModelError> {
    assert!(matches!(units, Units::Mom)); //Dejo un assert para que sepamos que esta mal. las units tienen
                                  //que volar.
                                  //
    let t = recursion_varying_anual(m,phi,dt).unwrap();
    let h = toeplitz_covariance_to_impulse(t);

    Ok(SigmaH {sigma: sigma, h: h})
}

//def Powerlaw(
//        *,
//        m,
//        kappa,
//        sigma,
//        units,
//        dt,
//        **kwargs
//    ):
//    d = -kappa/2.0
//    gmsv = gauss_markov_scale_variance(sigma=sigma,spectral_density=d,units=units,dt=dt)
//    return gmsv, recursion_Power_Flicker_RW(m,d)

//pub fn flicker() {
//}
//
//pub fn random_walk() {
//}
//
//pub fn ggm() {
//}
//
//pub fn varying_annual() {
//}
//
//pub fn matern() {
//}
//
//pub fn AR1() {
//}

fn recursion_power_flicker_rw(m: usize, d: f32) -> Result<Vec<f32>,ModelError> {
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
    Ok(h)
}

fn recursion_ggm(m: usize, d: f32, one_minus_phi: f32) -> Result<Vec<f32>,ModelError> {
    let mut h: Vec<f32> = vec![0.0; m];
    h[0] = 1.0;
    let mut h0: f32 = 1.0;
    
    for i in 1..m {
        h[i] = (d+i as f32-1.0)/i as f32 * h0 * (1.0 - one_minus_phi);
        h0 = (d+i as f32-1.0)/i as f32 * h0 * (1.0 - one_minus_phi);
    }
    Ok(h)
}

fn recursion_varying_anual(m: usize, phi: f32, dt: f32) -> Result<Vec<f32>,ModelError> {
    let mut t: Vec<f32> = vec![0.0; m];
    t[0] = 1.0/(2.0 * (1.0 - phi.powi(2)));
    let mut ti = t[0];
    for i in 1..m {
        let factor = (2.0 * PI * i as f32 * dt/365.25).cos();
        ti = ti * phi;
        t[i] = ti * factor;
    }
    Ok(t)
}

fn gauss_markov_scale_variance(
    sigma: f32, 
    spectral_density: f32, 
    units: Units, 
    dt: f32
    ) -> Result<f32,ModelError> {
    
    let sigma2: f32 = match units {
        Units::Mom => (dt as f32/365.25).powf(0.5*spectral_density),
        Units::Msf => (dt as f32/3600.0).powf(0.5*spectral_density),
    };
    Ok(sigma*sigma2)
}


#[cfg(test)]
mod tests {
    use crate::create_hs::{
        white, 
        recursion_varying_anual,
        recursion_power_flicker_rw,
        recursion_ggm,
        powerlaw,
        Units,
        SigmaH,
        varying_anual,
        flicker,
        random_walk,
        generalized_gauss_markov
    };

    #[test]
    fn test_white_ok() {
        let sigma = 1.0;
        let m: usize = 10;
        let sigma_h = white(m, sigma).unwrap();
        assert_eq!(sigma_h.sigma, sigma);
        assert_eq!(sigma_h.h, vec![sigma, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    }
    
    #[test]
    fn test_recursion_power_flicker_rw_ok() {
        // Write other tests
        let d: f32 = 0.5;
        let m: usize = 10;
        let response = recursion_power_flicker_rw(m, d).unwrap();
        // Taken from the python function (run the function in a shell)
        let expected: Vec<f32> = vec![1., 0.5, 0.375 , 0.3125, 0.2734375 ,0.24609375, 0.22558594, 0.20947266, 0.19638062, 0.18547058];
        assert_eq!(response, expected);
    }

    #[test]
    fn test_recursion_power_flicker_rw_ok2() {
        // Write other tests
        let d: f32 = 0.1;
        let m: usize = 100;
        let response = recursion_power_flicker_rw(m, d).unwrap();
        
        // Taken from the python function (run the function in a shell)
        let expected: Vec<f32> = vec!
            [1.        , 0.1       , 0.055     , 0.0385    , 0.0298375 ,
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

    #[test]
    fn test_recursive_ggn_ok(){
        let m: usize = 10;
        let d: f32 = 0.5;
        let one_minus_phi: f32 = 0.01;

        let response = recursion_ggm(m, d, one_minus_phi).unwrap();
        let expected: Vec<f32> = vec!
            [1., 0.495, 0.3675375 , 0.30321844, 0.26266297, 0.23403271, 
             0.21238468, 0.1952422 , 0.18120917, 0.16943057];
        
        assert_eq_float_vec!(response, expected, 0.000001);
    }

    //#[test]
    //fn gauss_markov_scale_variance_ok(){
    //    let sigma: f32 = 0.2;
    //    let spectral_density: f32 = 0.5;
    //    let units: Units = Units::Mom;
    //    
    //    let response = gauss_markov_scale_variance(0.2, 0.5, Units::Mom, 1);
    //    let expected: f32 = 0.04574908785331594;
    //    assert_eq_float!(response, expected, 0.001);
    //    
    //}
    //

    #[test]
    fn test_recursion_va_ok(){
        let m: usize = 30;
        let phi: f32 = 0.9;
        let dt: f32 = 45.0;

        let expected: Vec<f32> = vec!
          [ 2.63157895,  1.69352561,  0.04812309, -1.30981751, -1.72481892,
           -1.15902682, -0.0946563 ,  0.81698155,  1.12819177,  0.79031741,
            0.10336597, -0.50711688, -0.73642616, -0.53707432, -0.09475213,
            0.31307667,  0.47970378,  0.36382478,  0.07971134, -0.19210324,
           -0.31181833, -0.24573134, -0.06370291,  0.11705167,  0.20225416,
            0.16550532,  0.04919263, -0.07074449, -0.13089981, -0.11117542];

        let response = recursion_varying_anual(m, phi, dt).unwrap();
        println!("{:?}",response);
        println!("{:?}",expected);
        assert!(abs_diff_eq!(response[..], expected[..], epsilon= 0.0001));
    }

    #[test]
    fn test_va_ok(){
        let m: usize = 30;
        let phi: f32 = 0.9;
        let units = Units::Mom;
        let dt: f32 = 45.0;
        let sigma_in: f32 = 1.0;

        let expected_h: Vec<f32> = vec!
            [ 0.        ,  0.87786761,  0.76334962,  0.27141886, -0.26897556,
             -0.56604219, -0.51067108, -0.19878008,  0.15779781,  0.36411019,
              0.34082252,  0.14373641, -0.09106606, -0.23363574, -0.22694426,
             -0.10285035,  0.05144828,  0.14952681,  0.15077978,  0.07294884,
             -0.02824066, -0.09543655, -0.09995966, -0.05135322,  0.01487012,
              0.06074528,  0.06625521,  0.036664  , -0.00514098, -0.04780201];

        let SigmaH { sigma, h } = varying_anual(m, phi, units , dt, sigma_in).unwrap();

        println!("{:?}",h);
        println!("{:?}",expected_h);
        assert!(abs_diff_eq!(h[..], expected_h[..], epsilon= 0.0001));
        assert!(abs_diff_eq!(sigma, 1.0));
    }

    #[test]
    fn test_power_law_ok(){
        //create_hs.Powerlaw(m=20, kappa=-1, sigma=2, units="mom", dt=0.5)
        let m: usize = 20;
        let kappa: f32 = -1.0;
        let sigma: f32 = 2.0;
        let dt: f32 = 0.5;

        let expected_number: f32 = 0.3847024397698062;
        let expected_array: Vec<f32> = vec!
            [1. , 0.5, 0.375, 0.3125, 0.2734375 , 0.24609375, 0.22558594, 
             0.20947266, 0.19638062, 0.18547058, 0.17619705, 0.1681881 , 
             0.16118026, 0.15498102, 0.14944598, 0.14446445, 0.13994993, 
             0.13583376, 0.1320606 , 0.12858532];
        
        let response = powerlaw(m, kappa, sigma, Units::Mom, dt).unwrap();
        assert_eq!(response.sigma , expected_number);
        assert_eq_float_vec!(response.h, expected_array, 0.001);
    }

    #[test]
    // fn flicker(m: usize, sigma: f32, units: Units, dt: f32) -> Result<SigmaH, ModelError> {
    fn test_flicker_ok()
    {
        let m: usize = 30;
        let sigma: f32 = 1.0;
        let units = Units::Mom;
        let dt: f32 = 1.5;

        let expected_h: Vec<f32> = vec!
            [1.        , 0.5       , 0.375     , 0.3125    , 0.2734375 ,
             0.24609375, 0.22558594, 0.20947266, 0.19638062, 0.18547058,
             0.17619705, 0.1681881 , 0.16118026, 0.15498102, 0.14944598,
             0.14446445, 0.13994993, 0.13583376, 0.1320606 , 0.12858532,
             0.12537069, 0.12238567, 0.11960418, 0.11700409, 0.1145665 ,
             0.11227517, 0.11011603, 0.10807685, 0.10614691, 0.10431679];
        let expected_sigma: f32 = 0.2531484418502317;

        let SigmaH { sigma, h } = flicker(m, sigma, units , dt).unwrap();

        println!("{:?}",h);
        println!("{:?}",expected_h);

        assert!(abs_diff_eq!(h[..], expected_h[..], epsilon= 0.0001));
        assert!(abs_diff_eq!(sigma, expected_sigma));
    }

    #[test]
    //fn random_walk(m: usize, sigma: f32, units: Units, dt: f32) -> Result<SigmaH, ModelError> {
    fn test_random_walk_ok()
    {
        let m: usize = 30;
        let sigma: f32 = 1.0;
        let units = Units::Mom;
        let dt: f32 = 1.5;

        let expected_h: Vec<f32> = vec! [1. ;m];
        let expected_sigma: f32 = 0.06408413361120015;

        let SigmaH { sigma, h } = random_walk(m, sigma, units , dt).unwrap();

        println!("{:?}",h);
        println!("{:?}",expected_h);

        assert!(abs_diff_eq!(h[..], expected_h[..], epsilon= 0.0001));
        assert!(abs_diff_eq!(sigma, expected_sigma));
    }

    #[test]
    //fn generalized_gauss_markov(m: usize, kappa: f32, one_minus_phi: f32, sigma: f32, units: Units, dt: f32) -> Result<SigmaH, ModelError> {
    fn test_generalized_gauss_markov_ok()
    {
        let m: usize = 30;
        let kappa: f32 = 0.3;
        let one_minus_phi: f32 = -0.3;
        let sigma: f32 = 1.0;
        let units = Units::Mom;
        let dt: f32 = 1.5;

        let expected_h: Vec<f32> = vec!
            [ 1.        , -0.195     , -0.1077375 , -0.08636956, -0.07999981,
             -0.08007981, -0.08415053, -0.09142354, -0.10176583, -0.11539114,
             -0.13275751, -0.15454181, -0.18165102, -0.21525646, -0.25684708,
             -0.30830211, -0.37198577, -0.45086863, -0.54868208, -0.67011408,
             -0.82105728, -1.00892301, -1.243039  , -1.53515317, -1.90007187,
             -2.35646913, -2.9279129 , -3.64416714, -4.54284479, -5.67150675];

        let expected_sigma: f32 = 1.51003642143141;

        let SigmaH { sigma, h } = generalized_gauss_markov(
                                    m, kappa, one_minus_phi, sigma, units , dt
                                    ).unwrap();

        println!("{:?}",h);
        println!("{:?}",expected_h);

        assert!(abs_diff_eq!(h[..], expected_h[..], epsilon= 0.0001));
        assert!(abs_diff_eq!(sigma, expected_sigma));
    }

}
