use std::{collections::HashMap, env, fs};

use serde::{Deserialize, Serialize};
use serde_json::Value;
use strum::Display;

use crate::network::NetworkModel;

#[derive(Clone, Copy, Serialize, Deserialize, Display, Debug, clap::ValueEnum, PartialEq, Eq)]
pub enum OutputModel {
    AdjacencyList,
    AdjacencyMatrix,
    EdgeList,
    NetRustObject,
}

#[derive(Deserialize, Debug, Clone)]
#[serde(untagged)]
pub enum NetworkPars {
    BarabasiAlbert(BarabasiAlbertPars),
    Complete(CompletePars),
    Configuration(ConfigurationPars),
    ConfigurationCorrelated(ConfigurationPars),
    ConfigurationUncorrelated(ConfigurationPars),
    ErdosRenyi(ErdosRenyiPars),
    Lattice(LatticePars),
    LatticePBC(LatticePars),
    //RandomGeometric(RandomGeometricPars),
    Regular(RegularPars),
    WattsStrogatz(WattsStrogatzPars),
}

#[derive(Clone, Copy, Serialize, Deserialize, Debug, PartialEq)]
pub struct BarabasiAlbertPars {
    pub attachment_rate: Option<usize>,
    pub average_degree: Option<usize>,
    pub size: usize,
    pub size_initial_cluster: Option<usize>,
}

#[derive(Clone, Copy, Serialize, Deserialize, Debug, PartialEq)]
pub struct CompletePars {
    pub size: usize,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct ConfigurationPars {
    pub degree_minimum: usize,
    pub degree_maximum: usize,
    pub power_law_exponent: f64,
    pub size: usize,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct ErdosRenyiPars {
    pub average_degree: Option<usize>,
    pub probability_connection: Option<f64>,
    pub size: usize,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct LatticePars {
    pub nxcells: usize,
    pub nycells: usize,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct RandomGeometricPars {
    pub radius: f64,
    pub size: usize,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct RegularPars {
    pub average_degree: usize,
    pub size: usize,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct WattsStrogatzPars {
    pub average_degree: usize,
    pub probability_rewiring: f64,
    pub size: usize,
}

pub fn load_parameters(
    model_network: NetworkModel,
) -> Result<NetworkPars, Box<dyn std::error::Error>> {
    let current_dir = env::current_dir()?;
    let filename = "config_network.json";
    let path = current_dir.join("config").join(filename);

    let data = fs::read_to_string(path)?;
    let json: HashMap<String, Value> = serde_json::from_str(&data)?;

    let model_key = match model_network {
        NetworkModel::BarabasiAlbert => "BarabasiAlbert",
        NetworkModel::Complete => "Complete",
        NetworkModel::Configuration => "Configuration",
        NetworkModel::ConfigurationCorrelated => "ConfigurationCorrelated",
        NetworkModel::ConfigurationUncorrelated => "ConfigurationUncorrelated",
        NetworkModel::ErdosRenyi => "ErdosRenyi",
        NetworkModel::Lattice => "Lattice",
        NetworkModel::LatticePBC => "LatticePBC",
        //NetworkModel::RandomGeometric => "RandomGeometric",
        NetworkModel::Regular => "Regular",
        NetworkModel::WattsStrogatz => "WattsStrogatz",
    };

    let model_params = json.get(model_key).ok_or_else(|| {
        format!(
            "Model key {} not found in the configuration file",
            model_key
        )
    })?;

    let network_params: NetworkPars = serde_json::from_value(model_params.clone())?;

    Ok(network_params)
}
