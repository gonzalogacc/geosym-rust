use rand::Rng;
use rand_distr::{Normal};
use crate::fftconvolve::fftconvolve_vec;
use crate::create_hs::SigmaH;

// Crear el ruido a partir de un impulso (SigmaH) y un generador de numeros aleatorios.
// Pasar los generadores es una práctica que nos puede ayudar a futuro
// y en el presente es necesario para poder pasar generadores con semilla y testear.
fn simple_noise<R: Rng>(impulso: SigmaH, mut rng: R) -> Vec<f32> {
    let SigmaH { sigma, h } = impulso;
    let normal = Normal::new(0.0,1.0).unwrap();
    let m = h.len();

    let mut vals: Vec<f32> = (0..m).map(|_| rng.sample(&normal)).collect();

    for v in vals.iter_mut() {
        *v = *v * sigma;
    }

    fftconvolve_vec(h,vals)[0..m].to_vec()
}

//
// La función que ahora llamo simple_noise es un pedacito de adentro de la función create_noise de
// python.
// Create noise recibe una lista de modelos y devuelve la suma. En este caso podría ser
// tranquilamente una suma de llamadas a simple_noise.
//

#[cfg(test)]
mod tests {
    use rand_chacha::ChaCha8Rng;
    use rand::SeedableRng;
    use crate::create_hs::{powerlaw,white,Units};
    use crate::noise::simple_noise;

    use rand::Rng;
    use rand_distr::{Normal};
    #[test]
    fn test_random_numbers(){
        // Test para asegurarnos de que el PRNG no cambió.
        // Si este test falla van a fallar los que le siguen.
        
        let mut rng = ChaCha8Rng::seed_from_u64(42);
        let normal = Normal::new(0.0,1.0).unwrap();
        let vals: Vec<f32> = (0..20).map(|_| rng.sample(&normal)).collect();

        //Test for saved values.
        let expected_array: Vec<f32> = vec![0.47798123955726624, 1.3340705633163452, -0.21086668968200684, 0.4763469099998474, -0.5120906233787537, -0.9339784383773804, -1.0023778676986694, 0.9166635870933533, 2.1215765476226807, -0.7185473442077637, 0.031378861516714096, 1.0449801683425903, 2.032183885574341, 0.3553946912288666, 0.573213517665863, -0.8301297426223755, 0.2881774306297302, -0.5208272933959961, -0.024084562435746193, -0.8026301860809326];

      assert_eq!(vals, expected_array);

    }

    #[test]
    fn test_white_ok(){
        //create_hs.Powerlaw(m=20, kappa=-1, sigma=2, units="mom", dt=0.5)
        let m: usize = 20;
        let sigma: f32 = 5.4881350392732475;

        let expected_array: Vec<f32> = vec![ 2.623226, 7.321559, -1.157265, 2.614256, -2.810422, -5.125800, -5.501185, 5.030774, 11.643499, -3.943485, 0.172211, 5.734992, 11.152900, 1.950454, 3.145873, -4.555864, 1.581557, -2.858371, -0.132179, -4.404943];

        
        let impulso = white(m, sigma).unwrap();
        let rng = ChaCha8Rng::seed_from_u64(42);

        let response = simple_noise( impulso, rng);

        println!("{:?}",response);

        assert_eq_float_vec!(response, expected_array, 0.001);
    }

    #[test]
    fn test_power_law_ok(){
        //create_hs.Powerlaw(m=20, kappa=-1, sigma=2, units="mom", dt=0.5)
        let m: usize = 20;
        let kappa: f32 = 0.17953273198758746;
        let sigma: f32 = 1.0;
        let dt: f32 = 8.76127;

        let expected_array: Vec<f32> = vec![0.565097, 1.526487, -0.413966, 0.506407, -0.697515, -1.104593, -1.101335, 1.219118, 2.473125, -1.080844, 0.009084, 1.200685, 2.264048, 0.131247, 0.490321, -1.162136, 0.311531, -0.695711, -0.034496, -0.973195];
        
        let impulso = powerlaw(m, kappa, sigma, Units::Mom, dt).unwrap();
        let rng = ChaCha8Rng::seed_from_u64(42);

        let response = simple_noise( impulso, rng);


        assert_eq_float_vec!(response, expected_array, 0.001);
    }
}
