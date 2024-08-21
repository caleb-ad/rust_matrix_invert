use std::{cell::{RefCell, RefMut}, collections::{HashMap, VecDeque}, mem::MaybeUninit};
use num::Zero;
use std::rc::Rc;
use crate::matrix::{Matrix, SubMatrix};

#[derive(Debug)]
pub struct Graph<N, E> {
    nodes: Vec<N>,
    edges: E
}

pub trait Vertex {
    fn new(id: usize) -> Self;

    fn id(&self) -> usize;
}

pub trait Edge {
    fn new(size: usize) -> Self;

    fn connected(&self, a: usize, b: usize) -> bool;

    fn connect(&mut self, a: usize, b: usize);

    fn edges(&self, id: usize) -> impl Iterator<Item = usize>;
}

pub trait Visited{
    fn visited(&self) -> usize;

    fn visit(&mut self);
}

pub trait Distance {
    fn dist(&self) -> Option<usize>;

    fn set_dist(&mut self, d: usize);
}

pub trait Tree {
    fn parent(&self) -> Option<usize>;

    fn set_parent(&mut self, id: usize);
}

impl<N: Vertex, E: Edge> Graph<N, E> {
    pub fn new(size: usize) -> Self {
        Graph{nodes: (0..size).map(|i| Vertex::new(i)).collect(), edges: E::new(size)}
    }

    pub fn from_matrix<T: PartialEq + Zero>(m: SubMatrix<'_, T>) -> Self {
        assert!(m.dim.1.row == m.dim.1.col);
        let mut edges = E::new(m.dim.1.row);
        let vertices = (0..m.dim.1.row).map(|i| Vertex::new(i)).collect();

        for i in 0..m.dim.1.row { for j in 0..m.dim.1.col {
            if *m.get_ref((i, j)) != T::zero() {
                edges.connect(i, j)
            }
        }}

        Graph{nodes: vertices, edges}
    }
}

impl<N: Vertex + Visited, E: Edge> Graph<N, E> {
    fn unvisited(&self, i: usize) -> impl Iterator<Item = usize> + '_ {
        self.edges.edges(i).filter(|e| self.nodes[*e].visited() == 0)
    }

    fn visited(&self, i: usize) -> impl Iterator<Item = usize> + '_ {
        self.edges.edges(i).filter(|e| self.nodes[*e].visited() > 0)
    }
}

impl<N: Vertex + Tree, E: Edge> Graph<N, E> {
    fn parent(&self, i: usize) -> Option<&N> {
        self.nodes[i].parent().map(|n| &self.nodes[n])
    }

    fn not_parent(&self, i: usize) -> impl Iterator<Item = usize> + '_ {
        self.edges.edges(i).filter(move |v| *v != self.nodes[i].parent().unwrap_or(*v + 1))
    }
}

impl<N: Vertex + Visited + Distance + Tree, E: Edge> Graph<N, E> {
    // BFS to get cycle 'starts'
    // reverse history from each start until convergence
    pub fn get_cycles(&mut self) -> Vec<Vec<usize>> {
        let mut queue = VecDeque::new();
        let mut cycles = Vec::new();
        const START: usize = 0;
        self.nodes[START].set_dist(0);
        self.nodes[START].visit();
        queue.push_back(START);
        while let Some(pos) = queue.pop_front() {
            for n_id in self.not_parent(pos).collect::<Vec<usize>>() {
                // println!("next: {}", next);
                if self.nodes[n_id].visited() == 1 {
                    let mut cycle = vec![];
                    // println!("cycle at {} meets {}", pos, next);
                    self._collect_cycle(&mut cycle, pos, self.nodes[n_id].id());
                    cycles.push(cycle);
                } else if self.nodes[n_id].visited() == 0 {
                    let prev_dist = self.nodes[pos].dist().unwrap();
                    self.nodes[n_id].set_dist(prev_dist + 1);
                    self.nodes[n_id].visit();
                    self.nodes[n_id].set_parent(pos);
                    queue.push_back(n_id);
                    // self._get_cycles(queue, cycles, self.nodes[n_id].id());
                }
            }
            self.nodes[pos].visit()
        }
        cycles
    }

    fn _collect_cycle(&self, cycle: &mut Vec<usize>, pos1: usize, pos2: usize) {
        let d1 = self.nodes[pos1].dist().unwrap();
        let d2 = self.nodes[pos2].dist().unwrap();
        if pos1 == pos2 {
            cycle.push(pos1);
            return;
        }
        let (npos1, npos2) = if d1 > d2  {
            // println!("pos1 {:?} pos2 {:?}", pos1, pos2);
            cycle.push(pos1);
            (self.parent(pos1).unwrap().id(), pos2)
        } else {
            cycle.push(pos2);
            (pos1, self.parent(pos2).unwrap().id())
        };
        self._collect_cycle(cycle, npos1, npos2);
    }

    pub fn tree_partition(&mut self) -> Vec<Vec<usize>> {
        let cycles = self.get_cycles();
        let mut node_to_part: HashMap<usize, Rc<RefCell<Vec<usize>>>> = HashMap::new();
        for node in self.nodes.iter() {
            let r = Rc::new(RefCell::new(vec![node.id()]));
            node_to_part.insert(node.id(), r);
        }
        for cycle in cycles {
            // cycles should be sorted in distance
            let mut idx = 0;
            while idx < cycle.len() - 1 {
                if self.nodes[cycle[idx]].dist() == self.nodes[cycle[idx + 1]].dist() {
                    let part = node_to_part[&cycle[idx]].clone();
                    let mut part1 = part.borrow_mut();
                    match node_to_part[&cycle[idx + 1]].clone().try_borrow_mut() {
                        Ok(mut other_part) => {
                            part1.append(&mut other_part);
                            // update partition corresponding to node at idx + 1
                            node_to_part.insert(cycle[idx + 1], part.clone());
                        }
                        Err(_) => (), //nodes are already in the same partition
                    }
                    idx += 1;
                }
                idx += 1;
            }
        }
        node_to_part.into_values().filter_map(|rc| Rc::into_inner(rc).map(|r| r.into_inner())).collect()
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Cycle {
    idx: usize,
    prev: Option<usize>,
    visited: usize,
    dist: Option<usize>,
}

#[derive(Clone, Debug)]
pub struct Undirected {
    edges: Vec<Vec<bool>>,
}

impl Edge for Undirected {
    fn new(size: usize) -> Self {
        let mut t1 = Vec::new();
        let mut t2 = Vec::new();
        t2.resize(size, false);
        t1.resize(size, t2);
        Undirected{edges: t1}
    }

    fn connect(&mut self, a: usize, b: usize) {
        self.edges[a][b] = true;
        self.edges[b][a] = true;
    }

    fn connected(&self, a: usize, b: usize) -> bool {
        self.edges[a][b]
    }

    fn edges(&self, id: usize) -> impl Iterator<Item = usize> {
        self.edges[id].iter().enumerate().filter_map(|(idx, e)| if *e {Some(idx)} else {None})
    }
}

impl Vertex for Cycle {
    fn new(id: usize) -> Self{
        Cycle{idx: id, prev: None, visited: 0, dist: None}
    }

    fn id(&self) -> usize {
        self.idx
    }
}

impl Visited for Cycle {
    fn visit(&mut self) {
        self.visited += 1;
    }

    fn visited(&self) -> usize {
        self.visited
    }
}

impl Distance for Cycle {
    fn dist(&self) -> Option<usize> {
        return self.dist
    }

    fn set_dist(&mut self, d: usize) {
        self.dist = Some(d);
    }
}

impl Tree for Cycle {
    fn parent(&self) -> Option<usize> {
        self.prev
    }

    fn set_parent(&mut self, id: usize) {
        self.prev = Some(id);
    }
}


#[cfg(test)]
mod tests {
    use std::default;

    use super::*;

    #[test]
    fn test_get_cycles() {
        let mut g = Graph::<Cycle, Undirected>::new(14);
        g.edges.connect(0, 1);
        g.edges.connect(1, 2);
        g.edges.connect(0, 3);
        g.edges.connect(2, 0);
        g.edges.connect(2, 4);
        g.edges.connect(0, 5);
        g.edges.connect(5, 6);
        g.edges.connect(6, 7);
        g.edges.connect(7, 2);
        g.edges.connect(6, 8);
        g.edges.connect(8, 9);
        g.edges.connect(9, 10);
        g.edges.connect(10, 11);
        g.edges.connect(11, 12);
        g.edges.connect(12, 9);
        g.edges.connect(10, 13);
        g.edges.connect(11, 1);
        println!("{:?}", g.get_cycles());
    }

    #[test]
    fn test_tree_partition() {
        let mut g = Graph::<Cycle, Undirected>::new(14);
        g.edges.connect(0, 1);
        g.edges.connect(1, 2);
        g.edges.connect(0, 3);
        g.edges.connect(2, 0);
        g.edges.connect(2, 4);
        g.edges.connect(0, 5);
        g.edges.connect(5, 6);
        g.edges.connect(6, 7);
        g.edges.connect(7, 2);
        g.edges.connect(6, 8);
        g.edges.connect(8, 9);
        g.edges.connect(9, 10);
        g.edges.connect(10, 11);
        g.edges.connect(11, 12);
        g.edges.connect(12, 9);
        g.edges.connect(10, 13);
        g.edges.connect(11, 1);
        println!("{:?}", g.tree_partition());
    }

    #[test]
    fn test_from_matrix() {
        let a: Matrix<i32> = [[1,0,1],[0,1,0],[0,1,0]].into();
        let g: Graph<Cycle, Undirected> = Graph::from_matrix((&a).into());
        println!("{:?}", g);
        let a: Matrix<i32> = [[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,0,0,0]].into();
        let g: Graph<Cycle, Undirected> = Graph::from_matrix((&a).into());
        println!("{:?}", g);
    }

}
