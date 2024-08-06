use std::{collections::VecDeque, mem::MaybeUninit};

// Assuming Bidirectional
#[derive(Debug, Clone, Copy)]
pub struct Graph<C, const N: usize> {
    a: [[C; N]; N],
}

pub trait Vertex {
    fn id(&self) -> usize;
}

pub trait Edge {
    fn connected(&self) -> bool;

    fn connect(&mut self);
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
    fn parent(&self) -> bool;

    fn set_parent(&mut self, id: usize);
}


impl<C: Default + Copy, const N: usize> Default for Graph<C, N> {
    fn default() -> Self {
        Graph { a: [[C::default(); N]; N] }
    }
}

impl<C: Edge, const N: usize> Graph<C, N> {
    fn edges(&self, i: usize) -> impl Iterator<Item = &C> + '_ {
        self.a[i].iter().filter(| n| n.connected())
    }

    // create bidirectional edge
    fn connect(&mut self, v1: usize, v2: usize) {
        self.a[v1][v2].connect();
        self.a[v2][v1].connect();
    }

    fn new() -> Self {
        unsafe {
            let mut tmp: MaybeUninit<[[C;N];N]> = MaybeUninit::uninit();
            for i in 0..N {
                for j in 0..N {
                    (*tmp.as_mut_ptr())[i][j] = C::new(j);
                }
            }
            Graph{a: tmp.assume_init()}
        }
    }
}

impl<C: Edge + Visited, const N: usize> Graph<C, N> {
    fn unvisited(&self, i: usize) -> impl Iterator<Item = &C> + '_ {
        self.edges(i).filter(|e| !e.visited())
    }

    fn visited(&self, i: usize) -> impl Iterator<Item = &C> + '_ {
        self.edges(i).filter(|e| e.visited())
    }
}

impl<C: Tree + Edge, const N: usize> Graph<C, N> {
    fn make_parent(&mut self, parent: usize, child: usize) {
        self.a[child][parent].set_parent();
    }

    fn parent(&self, i: usize) -> Option<&C> {
        self.edges(i).find(|v| v.parent())
    }

    fn not_parent(&self, i: usize) -> impl Iterator<Item = &C> {
        self.edges(i).filter(|v| !v.parent())
    }
}


#[derive(Clone, Copy, Debug)]
pub struct Cycle {
    idx: usize,
    edge: bool,
    prev: bool,
    visited: Option<usize>,
    cycle: bool,
}

impl Edge for Cycle {
    fn new(idx: usize) -> Self {
        Cycle{idx: idx, edge: false, prev: false, visited: None, cycle: false}
    }

    fn idx(&self) -> usize {
        self.idx
    }

    fn connected(&self) -> bool {
        self.edge
    }

    fn connect(&mut self) {
        self.edge = true;
    }
}

impl Distance for Cycle {
    fn dist(&self) -> Option<usize> {
        return self.visited
    }

    fn set_dist(&mut self, d: usize) {
        self.visited = Some(d);
    }
}

impl Tree for Cycle {
    fn parent(&self) -> bool {
        self.prev
    }

    fn set_parent(&mut self) {
        self.prev = true;
    }
}

impl<const N: usize> Graph<Cycle, N> {
    // DFS to get cycle 'starts'
    // reverse history from each start until convergence
    fn get_cycles(&mut self) -> Vec<Vec<usize>> {
        let mut q = VecDeque::new();
        let mut c = Vec::new();
        const START: usize = 0;
        self.a[START][START].set_dist(0);
        self._get_cycles(&mut q, &mut c, START);
        c
    }

    fn _get_cycles(&mut self, queue: &mut VecDeque<usize>, cycles: &mut Vec<Vec<usize>>, pos: usize) {
        for next in self.not_parent(pos).map(|n| n.idx()).collect::<Vec<usize>>() {
            // println!("next: {}", next);
            if self.a[next][next].visited() && !self.a[next][next].cycle{
                let mut cycle = vec![];
                // println!("cycle at {} meets {}", pos, next);
                self._collect_cycle(&mut cycle, pos, next);
                cycles.push(cycle);
            } else if !self.a[next][next].visited() {
                self.a[next][next].set_dist(self.a[pos][pos].dist().unwrap() + 1);
                self.make_parent(pos, next);
                self._get_cycles(queue, cycles, next);
            }
        }
        self.a[pos][pos].cycle = true;
    }

    fn _collect_cycle(&self, cycle: &mut Vec<usize>, pos1: usize, pos2: usize) {
        let d1 = self.a[pos1][pos1].dist().unwrap();
        let d2 = self.a[pos2][pos2].dist().unwrap();
        if pos1 == pos2 {
            cycle.push(pos1);
            return;
        }
        let (npos1, npos2) = if d1 > d2  {
            // println!("pos1 {:?} pos2 {:?}", pos1, pos2);
            cycle.push(pos1);
            let npos1 = self.parent(pos1).unwrap().idx();
            (npos1, pos2)
        } else {
            cycle.push(pos2);
            let npos2 = self.parent(pos2).unwrap().idx();
            (pos1, npos2)
        };
        self._collect_cycle(cycle, npos1, npos2);
    }
}


#[cfg(test)]
mod tests {
    use std::default;

    use super::*;

    #[test]
    fn test_get_cycles() {
        let mut g = Graph::<Cycle, 14>::new();
        g.connect(0, 1);
        g.connect(1, 2);
        g.connect(0, 3);
        g.connect(2, 0);
        g.connect(2, 4);
        g.connect(0, 5);
        g.connect(5, 6);
        g.connect(6, 7);
        g.connect(7, 2);
        g.connect(6, 8);
        g.connect(8, 9);
        g.connect(9, 10);
        g.connect(10, 11);
        g.connect(11, 12);
        g.connect(12, 9);
        g.connect(10, 13);
        g.connect(11, 1);
        println!("{:?}", g.get_cycles());
    }

}