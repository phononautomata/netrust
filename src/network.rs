use rand::prelude::*;
use serde::{Deserialize, Serialize};
use serde_pickle::ser::SerOptions;
use std::{
    collections::{HashMap, HashSet, VecDeque},
    env,
    fs::{self, File},
    io::Write,
};
use strum::Display;
use uuid::Uuid;

use crate::{
    dists::build_powerlaw_degree_sequence,
    utils::{NetworkPars, OutputModel},
};

#[derive(Clone, Copy, Serialize, Deserialize, Display, Debug, clap::ValueEnum, PartialEq, Eq)]
pub enum NetworkModel {
    #[serde(rename = "complete")]
    Complete,
    #[serde(rename = "barabasi-albert")]
    BarabasiAlbert,
    #[serde(rename = "configuration")]
    Configuration,
    #[serde(rename = "configuration-correlated")]
    ConfigurationCorrelated,
    #[serde(rename = "configuration-uncorrelated")]
    ConfigurationUncorrelated,
    #[serde(rename = "erdos-renyi")]
    ErdosRenyi,
    #[serde(rename = "lattice")]
    Lattice,
    #[serde(rename = "lattice-pbc")]
    LatticePBC,
    //#[serde(rename = "random-geometric")]
    //RandomGeometric,
    #[serde(rename = "regular")]
    Regular,
    #[serde(rename = "watts-strogatz")]
    WattsStrogatz,
}

#[derive(Serialize)]
pub struct Node {
    pub id: usize,
    pub k: usize,
    pub neighbors: Vec<usize>,
}

impl Node {
    pub fn new(id: usize, k: usize, neighbors: Vec<usize>) -> Self {
        Self { id, k, neighbors }
    }
}

#[derive(Serialize)]
pub struct Network {
    nodes: Vec<Node>,
}

impl Default for Network {
    fn default() -> Self {
        Self::new()
    }
}

impl Network {
    pub fn new() -> Self {
        Self { nodes: Vec::new() }
    }

    pub fn inner(&self) -> &Vec<Node> {
        &self.nodes
    }

    pub fn inner_mut(&mut self) -> &mut Vec<Node> {
        &mut self.nodes
    }

    pub fn into_inner(self) -> Vec<Node> {
        self.nodes
    }

    pub fn add_node(&mut self, node: Node) {
        self.nodes.push(node);
    }

    pub fn add_edge(&mut self, id1: usize, id2: usize) {
        let node1_index = match self.nodes.iter().position(|node| node.id == id1) {
            Some(index) => index,
            None => return,
        };

        let node2_index = match self.nodes.iter().position(|node| node.id == id2) {
            Some(index) => index,
            None => return,
        };

        self.nodes[node1_index].neighbors.push(id2);
        self.nodes[node2_index].neighbors.push(id1);
    }

    pub fn average_degree(&self) -> f64 {
        let mut total_degree = 0;
        let nodes_count = self.nodes.len();

        for node in &self.nodes {
            total_degree += node.neighbors.len();
        }

        total_degree as f64 / nodes_count as f64
    }

    pub fn average_shortest_path_length(&self) -> f64 {
        let mut distance = vec![vec![std::usize::MAX; self.nodes.len()]; self.nodes.len()];
        for (i, node) in self.nodes.iter().enumerate() {
            distance[i][i] = 0;
            for &neighbor_id in &node.neighbors {
                distance[i][neighbor_id] = 1;
            }
        }
        for k in 0..self.nodes.len() {
            for i in 0..self.nodes.len() {
                for j in 0..self.nodes.len() {
                    distance[i][j] = distance[i][j].min(distance[i][k] + distance[k][j]);
                }
            }
        }
        let mut total_distance = 0;
        let mut count = 0;
        for (i, _) in self.nodes.iter().enumerate() {
            for (j, &dist) in distance[i].iter().enumerate() {
                if i != j && dist != std::usize::MAX {
                    total_distance += dist;
                    count += 1;
                }
            }
        }
        total_distance as f64 / count as f64
    }

    pub fn construct_string_network(&self, pars_net: &NetworkPars, uuid: &str) -> String {
        match pars_net {
            NetworkPars::BarabasiAlbert(pars) => {
                format!(
                    "ba_n{}_k{}_{}",
                    pars.size,
                    pars.average_degree.unwrap_or(pars.attachment_rate.unwrap()),
                    uuid
                )
            }
            NetworkPars::Complete(pars) => {
                format!("co_n{}", pars.size)
            }
            NetworkPars::ErdosRenyi(pars) => {
                format!(
                    "er_n{}_k{}_{}",
                    pars.size,
                    pars.average_degree.unwrap(),
                    uuid
                )
            }
            NetworkPars::Configuration(pars) => {
                format!(
                    "con_n{}_kmin{}_kmax{}_exp{}_{}",
                    pars.size,
                    pars.degree_minimum,
                    pars.degree_maximum,
                    pars.power_law_exponent,
                    uuid
                )
            }
            NetworkPars::ConfigurationCorrelated(pars) => {
                format!(
                    "cco_n{}_kmin{}_kmax{}_exp{}_{}",
                    pars.size,
                    pars.degree_minimum,
                    pars.degree_maximum,
                    pars.power_law_exponent,
                    uuid
                )
            }
            NetworkPars::ConfigurationUncorrelated(pars) => {
                format!(
                    "ucm_n{}_kmin{}_exp{}_{}",
                    pars.size, pars.degree_minimum, pars.power_law_exponent, uuid
                )
            }
            NetworkPars::Lattice(pars) => {
                format!("lat_nx{}_ny{}", pars.nxcells, pars.nycells)
            }
            NetworkPars::LatticePBC(pars) => {
                format!("lpb_nx{}_ny{}", pars.nxcells, pars.nycells)
            }
            NetworkPars::Regular(pars) => {
                format!("reg_n{}_k{}_{}", pars.size, pars.average_degree, uuid)
            }
            NetworkPars::WattsStrogatz(pars) => {
                format!(
                    "ws_n{}_k{}_p{}_{}",
                    pars.size, pars.average_degree, pars.probability_rewiring, uuid
                )
            }
        }
    }

    pub fn degree_sequence(&self) -> Vec<usize> {
        let mut degrees = vec![0; self.nodes.len()];
        for node in &self.nodes {
            degrees[node.id] = node.neighbors.len();
        }
        degrees
    }

    pub fn degree_distribution(&self) -> HashMap<usize, usize> {
        let mut distribution = HashMap::new();
        for node in &self.nodes {
            *distribution.entry(node.neighbors.len()).or_default() += 1;
        }
        distribution
    }

    pub fn maximum_degree(&self) -> usize {
        self.nodes
            .iter()
            .map(|node| node.neighbors.len())
            .max()
            .unwrap_or(0)
    }

    pub fn minimum_degree(&self) -> usize {
        self.nodes
            .iter()
            .map(|node| node.neighbors.len())
            .min()
            .unwrap_or(0)
    }

    pub fn moment_of_p_order(&self, p: usize) -> f64 {
        let mut total_degree_p = 0;
        let nodes_count = self.nodes.len();

        for node in &self.nodes {
            total_degree_p += node.neighbors.len().pow(p as u32);
        }

        total_degree_p as f64 / nodes_count as f64
    }

    pub fn number_of_nodes(&self) -> usize {
        self.nodes.len()
    }

    pub fn number_of_connections(&self) -> usize {
        let mut connections = 0;
        for node in &self.nodes {
            connections += node.neighbors.len();
        }
        connections
    }

    pub fn is_connected_bfs(&self) -> bool {
        let mut visited = vec![false; self.nodes.len()];
        let mut queue = VecDeque::new();
        queue.push_back(0);
        visited[0] = true;
        while let Some(node_id) = queue.pop_front() {
            for &neighbor_id in &self.nodes[node_id].neighbors {
                if !visited[neighbor_id] {
                    queue.push_back(neighbor_id);
                    visited[neighbor_id] = true;
                }
            }
        }
        visited.into_iter().all(|is_visited| is_visited)
    }

    pub fn is_connected_dfs(&self) -> bool {
        let mut visited = vec![false; self.nodes.len()];
        let mut stack = Vec::new();
        stack.push(0);
        visited[0] = true;
        while let Some(node_id) = stack.pop() {
            for &neighbor_id in &self.nodes[node_id].neighbors {
                if !visited[neighbor_id] {
                    stack.push(neighbor_id);
                    visited[neighbor_id] = true;
                }
            }
        }
        visited.into_iter().all(|is_visited| is_visited)
    }

    pub fn diameter(&self) -> usize {
        let mut distance = vec![vec![std::usize::MAX; self.nodes.len()]; self.nodes.len()];
        for (i, node) in self.nodes.iter().enumerate() {
            distance[i][i] = 0;
            for &neighbor_id in &node.neighbors {
                distance[i][neighbor_id] = 1;
            }
        }
        for k in 0..self.nodes.len() {
            for i in 0..self.nodes.len() {
                for j in 0..self.nodes.len() {
                    distance[i][j] = distance[i][j].min(distance[i][k] + distance[k][j]);
                }
            }
        }
    
        let mut diameter = 0;
        for (i, row) in distance.iter().enumerate().take(self.nodes.len()) {
            for (j, &dist) in row.iter().enumerate().take(self.nodes.len()) {
                if i != j {
                    diameter = diameter.max(dist);
                }
            }
        }
        diameter
    }

    pub fn clustering_coefficient(&self) -> f64 {
        let mut sum = 0.0;
        let mut count = 0;

        for i in 0..self.nodes.len() {
            let mut neighbors = HashSet::new();
            for &j in self.nodes[i].neighbors.iter() {
                neighbors.insert(j);
            }
            let k = self.nodes[i].neighbors.len();
            if k > 1 {
                let edges = k * (k - 1) / 2;
                let mut triangles = 0;
                for &j in self.nodes[i].neighbors.iter() {
                    for &k in self.nodes[j].neighbors.iter() {
                        if neighbors.contains(&k) {
                            triangles += 1;
                        }
                    }
                }
                triangles /= 2;
                sum += triangles as f64 / edges as f64;
                count += 1;
            }
        }

        if count == 0 {
            0.0
        } else {
            sum / count as f64
        }
    }

    pub fn generate_network(model_network: NetworkModel, pars_net: NetworkPars) -> Network {
        match model_network {
            NetworkModel::BarabasiAlbert => {
                if let NetworkPars::BarabasiAlbert(pars) = pars_net {
                    let average_degree = pars.average_degree;
                    let m = if let Some(m) = pars.attachment_rate {
                        m
                    } else {
                        (average_degree.unwrap() as f64 / 2.0) as usize
                    };
                    let m0 = if let Some(m0) = pars.size_initial_cluster {
                        m0
                    } else {
                        m
                    };

                    Network::generate_barabasi_albert(pars.size, m, m0)
                } else {
                    panic!("Expected NetworkPars::BarabasiAlbert");
                }
            }
            NetworkModel::Complete => {
                if let NetworkPars::Complete(pars) = pars_net {
                    Network::generate_complete_graph(pars.size)
                } else {
                    panic!("Expected NetworkPars::Complete");
                }
            }
            NetworkModel::Configuration => {
                if let NetworkPars::Configuration(pars) = pars_net {
                    let degrees = build_powerlaw_degree_sequence(
                        pars.size,
                        pars.power_law_exponent,
                        pars.degree_minimum,
                        pars.degree_maximum,
                    );
                    Network::generate_simple_configuration_model(degrees)
                } else {
                    panic!("Expected NetworkPars::Configuration");
                }
            }
            NetworkModel::ConfigurationCorrelated => {
                if let NetworkPars::ConfigurationCorrelated(pars) = pars_net {
                    let degrees = build_powerlaw_degree_sequence(
                        pars.size,
                        pars.power_law_exponent,
                        pars.degree_minimum,
                        pars.degree_maximum,
                    );
                    Network::generate_rewiring_configuration_model(degrees)
                } else {
                    panic!("Expected NetworkPars::ConfigurationCorrelated");
                }
            }
            NetworkModel::ConfigurationUncorrelated => {
                if let NetworkPars::ConfigurationUncorrelated(pars) = pars_net {
                    let degrees = build_powerlaw_degree_sequence(
                        pars.size,
                        pars.power_law_exponent,
                        pars.degree_minimum,
                        (pars.size as f64).sqrt() as usize,
                    );
                    Network::generate_rewiring_configuration_model(degrees)
                } else {
                    panic!("Expected NetworkPars::ConfigurationUncorrelated");
                }
            }
            NetworkModel::ErdosRenyi => {
                if let NetworkPars::ErdosRenyi(pars) = pars_net {
                    let average_degree = pars.average_degree;
                    let p = if let Some(connection_prob) = pars.probability_connection {
                        connection_prob
                    } else {
                        average_degree.unwrap() as f64 / (pars.size - 1) as f64
                    };

                    Network::generate_erdos_renyi(pars.size, p)
                } else {
                    panic!("Expected NetworkPars::ErdosRenyi");
                }
            }
            NetworkModel::Lattice => {
                if let NetworkPars::Lattice(pars) = pars_net {
                    Network::generate_lattice(pars.nxcells, pars.nycells)
                } else {
                    panic!("Expected NetworkPars::Lattice");
                }
            }
            NetworkModel::LatticePBC => {
                if let NetworkPars::LatticePBC(pars) = pars_net {
                    Network::generate_lattice_pbc(pars.nxcells, pars.nycells)
                } else {
                    panic!("Expected NetworkPars::LatticePBC");
                }
            }
            NetworkModel::Regular => {
                if let NetworkPars::Regular(pars) = pars_net {
                    Network::generate_k_regular_graph(pars.size, pars.average_degree)
                } else {
                    panic!("Expected NetworkPars::Regular");
                }
            }
            NetworkModel::WattsStrogatz => {
                if let NetworkPars::WattsStrogatz(pars) = pars_net {
                    Network::generate_watts_strogatz(
                        pars.size,
                        pars.average_degree,
                        pars.probability_rewiring,
                    )
                } else {
                    panic!("Expected NetworkPars::WattsStrogatz");
                }
            }
        }
    }

    pub fn generate_barabasi_albert(n: usize, m: usize, m0: usize) -> Network {
        if m > m0 {
            panic!("Number of edges per new node should be less than or total to initial number of nodes.");
        }
        if m >= n {
            panic!("Number of edges per new node must be less than total number of nodes.");
        }

        let mut nodes = Vec::new();
        //let mut edges = Vec::new();
        let mut cumulative_degrees = Vec::new();

        for i in 0..m0 {
            nodes.push(Node {
                id: i,
                k: 0,
                neighbors: vec![],
            });
            for j in 0..i {
                //edges.push((j, i));
                nodes[j].neighbors.push(i);
                nodes[i].neighbors.push(j);
                nodes[j].k += 1;
                nodes[i].k += 1;
            }
        }

        for node in &nodes {
            cumulative_degrees.push(node.k + cumulative_degrees.last().unwrap_or(&0));
        }

        let mut total_degree = cumulative_degrees.iter().sum();

        while nodes.len() < n {
            let next_node = nodes.len();
            nodes.push(Node {
                id: next_node,
                k: 0,
                neighbors: vec![],
            });

            let mut connected_nodes = Vec::new();
            let cumulative_current_degree = total_degree;

            while connected_nodes.len() < m {
                let rng = &mut rand::thread_rng();
                let target_value = rng.gen_range(0..=cumulative_current_degree);
                let target_node_idx = match cumulative_degrees.binary_search(&target_value) {
                    Ok(idx) => idx,
                    Err(idx) => idx,
                };

                if target_node_idx == next_node {
                    continue;
                }

                let target_node = nodes[target_node_idx].id;
                if !connected_nodes.contains(&target_node) {
                    connected_nodes.push(target_node);
                    //edges.push((target_node, next_node));
                    nodes[target_node].neighbors.push(next_node);
                    nodes[next_node].neighbors.push(target_node);
                    nodes[target_node].k += 1;
                    nodes[next_node].k += 1;
                    total_degree += 2;

                    for degree in cumulative_degrees.iter_mut().skip(target_node_idx) {
                        *degree += 1;
                    }
                }
            }

            cumulative_degrees.push(nodes[next_node].k + cumulative_degrees.last().unwrap_or(&0));
            total_degree = cumulative_degrees.iter().sum();
        }

        Network { nodes }
    }

    pub fn generate_complete_graph(n: usize) -> Self {
        let mut nodes = vec![];
        for i in 0..n {
            nodes.push(Node {
                id: i,
                k: 0,
                neighbors: vec![],
            });
        }

        for i in 0..n {
            for j in i + 1..n {
                nodes[i].neighbors.push(j);
                nodes[j].neighbors.push(i);
                nodes[i].k += 1;
                nodes[j].k += 1;
            }
        }

        Network { nodes }
    }

    pub fn generate_configuration_model(degrees: Vec<usize>) -> Self {
        let n = degrees.len();
        let mut nodes = vec![];
        let mut edges = vec![];

        for (i, &degree) in degrees.iter().enumerate().take(n) {
            nodes.push(Node {
                id: i,
                k: 0,
                neighbors: vec![],
            });
            for _ in 0..degree {
                edges.push(i);
            }
        }

        edges.shuffle(&mut rand::thread_rng());

        let mut next_edge_index = 0;

        for i in 0..n {
            let mut connected_neighbors = vec![];
            let mut j = degrees[i];

            while j > 0 {
                let neighbor = edges[next_edge_index];
                next_edge_index += 1;

                if neighbor == i || connected_neighbors.contains(&neighbor) {
                    continue;
                }

                nodes[i].neighbors.push(neighbor);
                nodes[i].k += 1;
                connected_neighbors.push(neighbor);
                j -= 1;
            }
        }

        Network { nodes }
    }

    pub fn generate_erdos_renyi(n: usize, p: f64) -> Self {
        let mut nodes = vec![];

        for i in 0..n {
            nodes.push(Node {
                id: i,
                k: 0,
                neighbors: vec![],
            });
        }

        for i in 0..n {
            for j in i + 1..n {
                if rand::random::<f64>() < p {
                    nodes[i].neighbors.push(j);
                    nodes[j].neighbors.push(i);
                    nodes[i].k += 1;
                    nodes[j].k += 1;
                }
            }
        }

        Network { nodes }
    }

    pub fn generate_k_regular_graph(n: usize, k: usize) -> Network {
        if n * k % 2 != 0 {
            panic!("Number of nodes * degree must be even.");
        }

        let mut nodes = Vec::new();
        let mut edges = Vec::new();
        let mut remaining_degrees = vec![k; n];

        for i in 0..n {
            for j in i + 1..n {
                if remaining_degrees[i] > 0 && remaining_degrees[j] > 0 {
                    edges.push((i, j));
                    remaining_degrees[i] -= 1;
                    remaining_degrees[j] -= 1;
                }
                if remaining_degrees[i] == 0 && remaining_degrees[j] == 0 {
                    break;
                }
            }
        }

        for i in 0..n {
            let mut k = 0;
            let mut neighbors = Vec::new();
            for &(x, y) in edges.iter() {
                if x == i {
                    neighbors.push(y);
                    k += 1;
                }
                if y == i {
                    neighbors.push(x);
                    k += 1;
                }
            }
            nodes.push(Node {
                id: i,
                k,
                neighbors,
            });
        }

        Network { nodes }
    }

    pub fn generate_lattice(nxcells: usize, nycells: usize) -> Network {
        let mut nodes = Vec::new();

        for y in 0..nycells {
            for x in 0..nxcells {
                let id = y * nxcells + x;
                let mut neighbors = Vec::new();

                if x > 0 {
                    neighbors.push(y * nxcells + (x - 1)); // left neighbor
                }
                if x < nxcells - 1 {
                    neighbors.push(y * nxcells + (x + 1)); // right neighbor
                }
                if y > 0 {
                    neighbors.push((y - 1) * nxcells + x); // top neighbor
                }
                if y < nycells - 1 {
                    neighbors.push((y + 1) * nxcells + x); // bottom neighbor
                }

                nodes.push(Node::new(id, neighbors.len(), neighbors));
            }
        }

        Network { nodes }
    }

    pub fn generate_lattice_pbc(nxcells: usize, nycells: usize) -> Network {
        let mut nodes = Vec::new();

        for y in 0..nycells {
            for x in 0..nxcells {
                let id = y * nxcells + x;
                let right = y * nxcells + (x + 1) % nxcells;
                let left = y * nxcells + (x + nxcells - 1) % nxcells;
                let up = ((y + 1) % nycells) * nxcells + x;
                let down = ((y + nycells - 1) % nycells) * nxcells + x;

                let neighbors = vec![right, left, up, down];
                nodes.push(Node::new(id, 4, neighbors));
            }
        }

        Network { nodes }
    }

    pub fn generate_simple_configuration_model(degrees: Vec<usize>) -> Self {
        let n = degrees.len();
        let mut nodes = vec![];
        let mut edges = vec![];

        for (i, &degree) in degrees.iter().enumerate().take(n) {
            nodes.push(Node {
                id: i,
                k: 0,
                neighbors: vec![],
            });
            for _ in 0..degree {
                edges.push(i);
            }
        }

        edges.shuffle(&mut rand::thread_rng());

        let mut next_edge_index = 0;
        for i in 0..n {
            for _ in 0..degrees[i] {
                let j = edges[next_edge_index];
                nodes[i].neighbors.push(j);
                nodes[i].k += 1;
                next_edge_index += 1;
            }
        }

        Network { nodes }
    }

    pub fn generate_rewiring_configuration_model(degrees: Vec<usize>) -> Self {
        let n = degrees.len();
        let mut nodes = vec![];
        let mut edges = vec![];

        for (i, &degree) in degrees.iter().enumerate().take(n) {
            nodes.push(Node {
                id: i,
                k: 0,
                neighbors: vec![],
            });
            for _ in 0..degree {
                edges.push(i);
            }
        }

        edges.shuffle(&mut rand::thread_rng());

        let mut bad_connections = vec![];
        let mut i = 0;
        while i <= edges.len() - 2 {
            if edges[i] == edges[i + 1] || nodes[edges[i]].neighbors.contains(&edges[i + 1]) {
                bad_connections.push((edges[i], edges[i + 1]));
            } else {
                nodes[edges[i]].neighbors.push(edges[i + 1]);
                nodes[edges[i + 1]].neighbors.push(edges[i]);
                nodes[edges[i]].k += 1;
                nodes[edges[i + 1]].k += 1;
            }
            i += 2;
        }

        let mut network = Network { nodes };

        network.swap_assitance(bad_connections);

        network
    }

    pub fn generate_watts_strogatz(n: usize, k: usize, p: f64) -> Network {
        let mut nodes = Vec::new();
        let mut edges = Vec::new();
        let mut degrees = vec![0; n];

        for i in 0..n {
            nodes.push(Node {
                id: i,
                k,
                neighbors: vec![],
            });
            for j in i + 1..i + k + 1 {
                let idx = j % n;
                edges.push((i, idx));
                nodes[i].neighbors.push(idx);
                nodes[idx].neighbors.push(i);
                nodes[i].k += 1;
                nodes[idx].k += 1;
                degrees[i] += 1;
                degrees[idx] += 1;
            }
        }

        for i in 0..n {
            for j in 0..k {
                let idx = (i + j + 1) % n;
                if rand::thread_rng().gen::<f64>() < p {
                    let mut neighbor = (rand::thread_rng().gen::<usize>() % n) as i32; // why i32???

                    while neighbor == i as i32 || nodes[i].neighbors.contains(&(&i) as &usize) {
                        neighbor = (rand::thread_rng().gen::<usize>() % n) as i32;
                    }

                    edges.retain(|edge| edge != &(i, idx) && edge != &(idx, i));
                    nodes[i].neighbors.retain(|neighbor| *neighbor != idx);
                    nodes[idx].neighbors.retain(|neighbor| *neighbor != i);
                    nodes[i].k -= 1;
                    nodes[idx].k -= 1;
                    degrees[i] -= 1;
                    degrees[idx] -= 1;
                    edges.push((i, neighbor as usize));
                    edges.push((neighbor as usize, i));
                    nodes[i].neighbors.push(neighbor as usize);
                    nodes[neighbor as usize].neighbors.push(i);
                    nodes[i].k += 1;
                    nodes[neighbor as usize].k += 1;
                    degrees[i] += 1;
                    degrees[neighbor as usize] += 1;
                }
            }
        }

        Network { nodes }
    }

    pub fn giant_component_size(&self) -> usize {
        let mut visited = vec![false; self.nodes.len()];
        let mut queue = VecDeque::new();
        let mut sizes = vec![];
        for node_id in 0..self.nodes.len() {
            if visited[node_id] {
                continue;
            }
            let mut size = 0;
            queue.push_back(node_id);
            visited[node_id] = true;
            while let Some(node_id) = queue.pop_front() {
                size += 1;
                for &neighbor_id in &self.nodes[node_id].neighbors {
                    if !visited[neighbor_id] {
                        queue.push_back(neighbor_id);
                        visited[neighbor_id] = true;
                    }
                }
            }
            sizes.push(size);
        }
        *sizes.iter().max().unwrap_or(&0)
    }

    pub fn repeated_connections(&self) -> usize {
        let mut count = 0;
        for node in &self.nodes {
            let mut seen = HashSet::new();
            for neighbor_id in &node.neighbors {
                if seen.contains(neighbor_id) {
                    count += 1;
                } else {
                    seen.insert(neighbor_id);
                }
            }
        }
        count / 2
    }

    pub fn save_to_file(&self, pars_net: &NetworkPars, model_output: OutputModel) {
        let uuid = Uuid::new_v4().to_string();
        let string_network = self.construct_string_network(pars_net, &uuid);

        match model_output {
            OutputModel::AdjacencyList => {
                let adj_list = self.to_adjacency_list();
                self.to_json(&adj_list, &string_network)
            }
            OutputModel::AdjacencyMatrix => {
                let adj_matrix = self.to_adjacency_matrix();
                self.to_json(&adj_matrix, &string_network)
            }
            OutputModel::EdgeList => {
                let edge_list = self.to_edge_list();
                self.to_json(&edge_list, &string_network)
            }
            OutputModel::NetRustObject => self.to_pickle(&string_network),
        }
    }

    pub fn second_moment(&self) -> f64 {
        let mut total_degree_squared = 0;
        let nodes_count = self.nodes.len();

        for node in &self.nodes {
            total_degree_squared += node.neighbors.len() * node.neighbors.len();
        }

        total_degree_squared as f64 / nodes_count as f64
    }

    pub fn self_connections(&self) -> usize {
        let mut count = 0;
        for node in &self.nodes {
            count += node
                .neighbors
                .iter()
                .filter(|&neighbor_id| *neighbor_id == node.id)
                .count();
        }
        count
    }

    pub fn size(&self) -> usize {
        self.nodes.len()
    }

    pub fn swap_assitance(&mut self, bad_connections: Vec<(usize, usize)>) {
        let mut rng = rand::thread_rng();
        let mut i = 0;
        let mut l = 0;
        let mut m = 0;

        while i < (bad_connections.len() - 1) {
            let mut search_candidates = true;
            while search_candidates {
                let mut empty_list = true;
                let mut number_of_neighs = 0;
                while empty_list {
                    l = rng.gen_range(0..self.size()) as usize;
                    number_of_neighs = self.nodes[l].neighbors.len();
                    if number_of_neighs != 0 {
                        empty_list = false;
                    }
                }
                let idx = rng.gen_range(0..number_of_neighs) as usize;
                m = self.nodes[l].neighbors[idx];

                let id1 = bad_connections[i].0 == l;
                let id2 = bad_connections[i].1 == m;
                let ex1 = self.nodes[l].neighbors.contains(&bad_connections[i].0);
                let ex2 = self.nodes[m].neighbors.contains(&bad_connections[i].1);

                if !id1 && !id2 && !ex1 && !ex2 {
                    search_candidates = false;
                }
            }

            // Now candidates fulfill conditions: swap connections
            // Remove the connection between m and l
            self.nodes[l].neighbors.retain(|&x| x != m);
            self.nodes[m].neighbors.retain(|&x| x != l);
            // Make a connection between bad_connection[i].0 and l and viceversa
            self.nodes[bad_connections[i].0].neighbors.push(l);
            self.nodes[l].neighbors.push(bad_connections[i].0);
            // Make a connection between bad_connection[i].1 and m and viceversa
            self.nodes[bad_connections[i].1].neighbors.push(m);
            self.nodes[m].neighbors.push(bad_connections[i].1);

            i += 1;
        }
    }

    pub fn to_adjacency_list(&self) -> HashMap<usize, Vec<usize>> {
        let mut adj_list = HashMap::new();
        for node in &self.nodes {
            adj_list.insert(node.id, node.neighbors.clone());
        }
        adj_list
    }

    pub fn to_adjacency_matrix(&self) -> Vec<Vec<bool>> {
        let mut matrix = vec![vec![false; self.nodes.len()]; self.nodes.len()];
        for node in &self.nodes {
            for &neighbor in &node.neighbors {
                matrix[node.id][neighbor] = true;
            }
        }
        matrix
    }

    pub fn to_edge_list(&self) -> Vec<(usize, usize)> {
        let mut edge_list = Vec::new();
        for node in &self.nodes {
            for &neighbor in &node.neighbors {
                edge_list.push((node.id, neighbor));
            }
        }
        edge_list
    }

    pub fn to_json<T: Serialize>(&self, data: &T, filename: &str) {
        let json = serde_json::to_string_pretty(data).expect("Failed to serialize data");
        let mut file = File::create(filename).expect("Failed to create file");
        file.write_all(json.as_bytes())
            .expect("Failed to write to file");
    }

    pub fn to_pickle(&self, string_network: &str) {
        let mut path = env::current_dir()
            .expect("Failed to get current directory")
            .join("data")
            .join("networks");

        if !path.exists() {
            fs::create_dir_all(&path).expect("Failed to create directory");
        }

        path.push(format!("net_nro_{}.pickle", string_network));

        let serialized = serde_pickle::to_vec(self, SerOptions::new()).unwrap();
        std::fs::write(path, serialized).unwrap();
    }
}

pub struct SpatialNode {
    pub id: usize,
    pub x: f64,
    pub y: f64,
    pub neighbors: Vec<usize>,
}

pub struct SpatialNetwork {
    nodes: Vec<SpatialNode>,
}

impl Default for SpatialNetwork {
     fn default() -> Self {
         Self::new()
     }
 }

impl SpatialNetwork {
    pub fn new() -> Self {
        Self { nodes: Vec::new() }
    }

    pub fn inner(&self) -> &Vec<SpatialNode> {
        &self.nodes
    }

    pub fn inner_mut(&mut self) -> &mut Vec<SpatialNode> {
        &mut self.nodes
    }

    pub fn into_inner(self) -> Vec<SpatialNode> {
        self.nodes
    }

    pub fn size(&self) -> u32 {
        self.nodes.len() as u32
    }

    pub fn add_node(&mut self, node: SpatialNode) {
        self.nodes.push(node);
    }

    pub fn add_edge(&mut self, id1: usize, id2: usize) {
        let node1_index = match self.nodes.iter().position(|node| node.id == id1) {
            Some(index) => index,
            None => return,
        };

        let node2_index = match self.nodes.iter().position(|node| node.id == id2) {
            Some(index) => index,
            None => return,
        };

        self.nodes[node1_index].neighbors.push(id2);
        self.nodes[node2_index].neighbors.push(id1);
    }

    pub fn build_random_geometric_network(n: usize, radius: f64) -> Self {
        let mut nodes = vec![];
        let mut rng = thread_rng();

        for i in 0..n {
            let id = i;
            let x = rng.gen_range(0.0..1.0);
            let y = rng.gen_range(0.0..1.0);
            nodes.push(SpatialNode {
                id,
                x,
                y,
                neighbors: vec![],
            });
        }

        for i in 0..n {
            for j in i + 1..n {
                let distance =
                    ((nodes[i].x - nodes[j].x).powi(2) + (nodes[i].y - nodes[j].y).powi(2)).sqrt();

                if distance <= radius {
                    nodes[i].neighbors.push(j);
                    nodes[j].neighbors.push(i);
                }
            }
        }

        SpatialNetwork { nodes }
    }
}
