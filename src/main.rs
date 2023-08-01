mod utils;
#[macro_use]
mod macros;
//mod simulatenoise;
mod create_hs;
mod onenoise;
mod fftconvolve;
use clap::Parser;
use std::error::Error;

#[derive(Parser)]
struct Cli {
    control_file: std::path::PathBuf,
}

fn main() -> Result<(), Box<dyn Error>> {
    
    println!("***************************************");
    println!("    simulatenoise, version 0.1.1");
    println!("***************************************");

    let args = Cli::parse();
    println!("{:?}", args.control_file);
    
 //   let config = utils::load_config_toml(args.control_file.as_path())?;
 //   println!("{:#?}", config);
    
 //   //let noise = simulatenoise::simulate_noise(&config);
 //   let wnoise: SigmaH = white_noise(10, 0.5).unwrap();
 //   println!("{:?}", wnoise.h);

 //   //for model in config.models(){
 //   //    println!("XXXXXXXXXXXXXXXXXXXXXXXXXXXXx");
 //   //    println!("{:#?}", model);

 //   //    let jj = config.get_model(&model).unwrap();
 //   //    println!("{:#?}", jj);

 //   //}
 //   //

 //   //println!("{:#?}", config.get_model("GGM")


 //   //let config = match utils::load_config_toml(args.control_file.as_path()) {
 //   //    Ok(c) => c,
 //   //    Err(e) => return Err("Error cargando archivo de control"),
 //   //};

 //   //let noise = simulatenoise::simulate_noise(&config);
 //   //println!("{:?}", noise);
 //   
    Ok(())
}
//
//    simulaciones = simulate_noise(control)
//
//    for simulacion in simulaciones:
//        print( '\n'.join(
//            map('{:.3f}'.format,
//                simulacion)
//                        )
//              )
