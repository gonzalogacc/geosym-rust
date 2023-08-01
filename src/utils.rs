use std::fs;
use serde::Deserialize;
use std::path::Path;
use std::io;
use crate::create_hs::{white, Units};


#[allow(dead_code)]
#[derive(Deserialize, Debug)]
pub struct Config {
    pub general: General,
    pub white_noise: WhiteNoise,
    pub ggm: Ggm,
}

//TODO: Add toher models configuration
#[allow(dead_code)]
#[derive(Deserialize, Debug)]
pub struct General {
    pub number_of_points: i32,
    pub sampling_period: f32,
    pub time_noise_start: i32,
    pub number_of_simulations: i32,
    pub repeatable_noise: bool,
    pub deterministic_noise: bool,
    pub order: Vec<String>
}

#[allow(dead_code)]
#[derive(Deserialize, Debug, Copy, Clone)]
pub struct Ggm {
    pub m: usize,
    pub sigma: f32,
    pub kappa: f32,
    pub one_minus_phi: f32,
    pub dt: f32, 
    pub units: Units
}

#[allow(dead_code)]
#[derive(Deserialize, Debug, Copy, Clone)]
pub struct WhiteNoise {
    m: usize,
    sigma: f32,
}

#[derive(Debug)]
enum Model {
    GGM(Ggm),
    WHITENOISE(WhiteNoise),
}

impl Model {
    fn call(&self) -> Vec<f32> {
        // Take the self type and match it to the correct model function
        let response = match &self {
            Model::WHITENOISE(x) => white(x.m, x.sigma).unwrap().h,
            Model::GGM(x) => white(x.m, x.sigma).unwrap().h,
        };
        response
    }
}

fn model_list_from_config(config: &Config) -> Vec<Model> {
    let mut res: Vec<Model> = Vec::new();
    for model in config.general.order.iter() {
        match model.as_str() {
            "GGM" => res.push(Model::GGM(config.ggm)),
            "WHITENOISE" => res.push(Model::WHITENOISE(config.white_noise)),
            _ => panic!("No model name")
        }
    }
    res
}

pub fn load_config_toml(filename: &Path) -> Result<Config, io::Error> {
    let content = fs::read_to_string(filename)?; 
    let config: Config = toml::from_str(&content).expect("Error parsing file");
    let lista = model_list_from_config(&config);
    println!("---> {:#?}", lista);
    Ok(config)
}


#[cfg(test)]
mod tests {
    use crate::utils::{load_config_toml, model_list_from_config};
    // use crate::create_hs::Units;
    use std::path::Path;
    
    #[test]
    fn test_utils_get_function_list_ok() {
        let path = Path::new("./test2.toml");
        let config = load_config_toml(&path).unwrap();
        println!("{:#?}", config);
        
        //convert the config to a list of Models to be called
        let lista = model_list_from_config(&config);
        println!("---> {:#?}", lista);
        
        // call each model
        let mut response: Vec<Vec<f32>> = Vec::new();
        for i in lista.iter() {
            response.push(i.call())
        }

        println!("{:?}", response);
        //let expected: Vec<Vec<f32>> = vec![vec![1.0, 5.0], vec![1.0, 5.0]];
        //assert_eq!(response, expected);
        //assert_eq!(1,2);
    }

    //#[test]
    //fn test_utils_load_config_toml_ok() {
    //    let path = Path::new("./test2.toml");
    //    let config = load_config_toml(&path).unwrap();
    //    println!("{:#?}", config);

    //    match config.ggm {
    //        Some(x) => assert_eq!(x.m, 1000),
    //        None => panic!("Missing ggm"),
    //    }

    //    match config.linear {
    //        Some(x) => assert_eq!(x.m, 1000),
    //        None => panic!("Missing linear"),
    //    }
    //}
}
