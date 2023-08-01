use std::io;
//use std::error::Error;
//use crate::utils::Units;
use serde::Deserialize;

#[allow(dead_code)] // Puede que haya errores listados no implementados
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

type ImpulseH = Vec<f32>;

#[derive(Debug)]
pub struct SigmaH {
    pub sigma: f32,
    pub h: ImpulseH,
}

impl SigmaH {
    fn new(sigma: f32, h: ImpulseH) -> SigmaH {
        SigmaH {
            sigma,
            h,
        }
    }
}

pub fn white(m: usize, sigma: f32) -> Result<SigmaH, ModelError> {
    // Create impulse function for White noise.
    let mut h: ImpulseH = vec![0.0; m];
    h[0] = 1.0; //sigma; El sigma se aplica en create_noise
    Ok(SigmaH::new(sigma, h))
}

pub fn powerlaw(m: usize, kappa: f32, sigma: f32, units: Units, dt: f32) -> Result<SigmaH, ModelError> {
    let d = -kappa/2.0;
    let gmsv = gauss_markov_scale_variance(sigma, d, units, dt).unwrap();
    let rpf = recursion_power_flicker_rw(m, d).unwrap();
    Ok(SigmaH::new(gmsv,rpf))
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

fn recursion_power_flicker_rw(m: usize, d: f32) -> Result<ImpulseH, ModelError> {
    // Recursion to create impulse function for Powerlay, Flicker or RW noise
    // Flicker is Powerlaw with spectral density 0.5
    // RandomWalk is Powerlaw with spectral density 1.0
    let mut h: ImpulseH = vec![0.0; m];
    
    h[0] = 1.0;
    let mut h0: f32 = 1.0;
    for i in 1..m {
        h[i] = (d+i as f32-1.0)/i as f32 * h0;
        h0 = (d+i as f32-1.0)/i as f32 * h0;
    }
    Ok(h)
}

fn recursion_ggm(m: usize, d: f32, one_minus_phi: f32) -> Result<ImpulseH, ModelError> {
    let mut h: Vec<f32> = vec![0.0; m];
    h[0] = 1.0;
    let mut h0: f32 = 1.0;
    
    for i in 1..m {
        h[i] = (d+i as f32-1.0)/i as f32 * h0 * (1.0 - one_minus_phi);
        h0 = (d+i as f32-1.0)/i as f32 * h0 * (1.0 - one_minus_phi);
    }
    Ok(h)
}

fn gauss_markov_scale_variance(
    sigma: f32, 
    spectral_density: f32, 
    units: Units, 
    dt: f32
    ) -> Result<f32, ModelError> {
    
    let sigma2: f32 = match units {
        Units::Mom => (dt as f32/365.25).powf(0.5*spectral_density),
        Units::Msf => (dt as f32/3600.0).powf(0.5*spectral_density),
    };
    Ok(sigma*sigma2)
}

//def gauss_markov_scale_variance(
//        *,
//        sigma,
//        spectral_density,
//        units,
//        dt,
//    ):
//    """
//    Gauss Markov Models needs scaling of the variance for taking into account
//    the time and period units.
//
//    This scale applies to Powerlaw noise, Flicker Noise, Random Walk Noise and
//    Generalized Gauss Markov noise.
//    """
//
//    if units=='mom':
//        sigma *= math.pow(dt/365.25,0.5*spectral_density)  # Already adjust for scaling
//    elif units=='msf':
//        sigma *= math.pow(dt/3600.0,0.5*spectral_density)  # Already adjust for scaling
//    else:
//        raise ValueError('unknown scaling: {0:s}'.format(units))
//
//    return sigma


#[cfg(test)]
mod tests {
    use crate::create_hs::{
        white, 
        recursion_power_flicker_rw,
        recursion_ggm,
        powerlaw,
        Units
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

    #[test]
    fn test_recursive_ggn_ok(){
        let m: usize = 10;
        let d: f32 = 0.5;
        let one_minus_phi: f32 = 0.01;

        let response = recursion_ggm(m, d, one_minus_phi).unwrap();
        let expected: Vec<f32> = vec![1., 0.495, 0.3675375 , 0.30321844, 0.26266297, 0.23403271, 0.21238468, 0.1952422 , 0.18120917, 0.16943057];
        
        assert_eq_float_vec!(response, expected, 0.000001);
    }

    //#[test]
    //fn gauss_markov_scale_variance_ok(){
    //    let sigma: f32 = 0.2;
    //    let spectral_density: f32 = 0.5;
    //    let units: Units = Units::mom;
    //    
    //    let response = gauss_markov_scale_variance(0.2, 0.5, Units::mom, 1);
    //    let expected: f32 = 0.04574908785331594;
    //    assert_eq_float!(response, expected, 0.001);
    //    
    //}

    #[test]
    fn test_power_law_ok(){
        //create_hs.Powerlaw(m=20, kappa=-1, sigma=2, units="mom", dt=0.5)
        let m: usize = 20;
        let kappa: f32 = -1.0;
        let sigma: f32 = 2.0;
        let dt: f32 = 0.5;

        let expected_number: f32 = 0.3847024397698062;
        let expected_array: Vec<f32> = vec![1. , 0.5, 0.375, 0.3125, 0.2734375 , 0.24609375, 0.22558594, 0.20947266, 0.19638062, 0.18547058, 0.17619705, 0.1681881 , 0.16118026, 0.15498102, 0.14944598, 0.14446445, 0.13994993, 0.13583376, 0.1320606 , 0.12858532];
        
        let response = powerlaw(m, kappa, sigma, Units::Mom, dt).unwrap();
        assert_eq!(response.sigma, expected_number);
        assert_eq_float_vec!(response.h, expected_array, 0.001);
    }

}


