use std::ops::Not;


#[derive(Clone, Copy, Debug, PartialEq)]
enum Dir {
    Left,
    Right,
    NoneLeft,
    NoneRight
}

impl From<Dir> for isize {
    fn from(value: Dir) -> Self {
        match value {
            Dir::Left => -1,
            Dir::Right => 1,
            Dir::NoneLeft | Dir::NoneRight => 0,
        }
    }
}

impl Not for Dir {
    type Output = Self;
    fn not(self) -> Self::Output {
        match self {
            Dir::Left => Dir::NoneLeft,
            Dir::NoneLeft => Dir::Right,
            Dir::Right => Dir::NoneRight,
            Dir::NoneRight => Dir::Left
        }
    }
}

//TODO: these can all be boxed slices
pub struct Permuter {
    p: Vec<usize>,
    val_to_idx: Vec<usize>,
    dir: Vec<Dir>,
    to_move: usize
}

impl Permuter {
    pub fn new(size: usize) -> Self {
        let mut tmp = Vec::new();
        tmp.resize(size, Dir::NoneRight);
        tmp[0] = Dir::Left;
        let mut tmp1: Vec<usize> = (0..size).collect();
        tmp1.swap(0, 1);

        Permuter{p: tmp1.clone(), val_to_idx: tmp1, dir: tmp, to_move: 0}
    }

    pub fn next(&mut self) -> Option<&[usize]> {
        while let Dir::NoneLeft | Dir::NoneRight = self.dir[self.to_move] {
            self.to_move = self.to_move.checked_sub(1)?;
        }

        // println!("{:?} {:?} {:?} {:?}", self.p, self.val_to_idx, self.dir, self.to_move);
        let old_pos = self.val_to_idx[self.to_move];
        let new_pos = (self.val_to_idx[self.to_move] as isize + Into::<isize>::into(self.dir[self.to_move])) as usize;
        self.val_to_idx.swap(self.to_move, self.p[new_pos]);
        self.p.swap(old_pos, new_pos);
        if new_pos == 0 || new_pos == self.p.len() - 1  {
            self.dir[self.to_move] = !self.dir[self.to_move];
        }
        if (old_pos == 0 && self.dir[self.p[old_pos]] == Dir::Left) || (old_pos == self.p.len() - 1 && self.dir[self.p[old_pos]] == Dir::Right) {
            self.dir[self.p[old_pos]] = !self.dir[self.p[old_pos]];
        }

        for i in self.to_move + 1..self.p.len() {
            self.dir[i] = !self.dir[i];
        }

        self.to_move = self.p.len() - 1;
        Some(&self.p)
    }
}

#[cfg(test)]
mod tests {
    use std::default;

    use super::*;

    #[test]
    fn test_permute() {
        let permutes = [
            &[0,1,2,3],
            &[0,1,3,2],
            &[0,3,1,2],
            &[3,0,1,2],
            &[3,0,2,1],
            &[0,3,2,1],
            &[0,2,3,1],
            &[0,2,1,3],
            &[2,0,1,3],
            &[2,0,3,1],
            &[2,3,0,1],
            &[3,2,0,1],
            &[3,2,1,0],
            &[2,3,1,0],
            &[2,1,3,0],
            &[2,1,0,3],
            &[1,2,0,3],
            &[1,2,3,0],
            &[1,3,2,0],
            &[3,1,2,0],
            &[3,1,0,2],
            &[1,3,0,2],
            &[1,0,3,2],
            &[1,0,2,3]
            ];
        let mut p = Permuter::new(4);
        for permutation in permutes {
            assert_eq!(permutation, p.next().unwrap());
        }
        assert!(p.next().is_none());
    }

    #[test]
    fn test_amnt() {
        let mut p = Permuter::new(5);
        for _ in 0..120 {assert!(p.next().is_some());}
        assert!(p.next().is_none());
        let mut p = Permuter::new(6);
        for _ in 0..720 {assert!(p.next().is_some());}
        assert!(p.next().is_none());
        let mut p = Permuter::new(12);
        for _ in 0..479001600_usize {assert!(p.next().is_some());}
        assert!(p.next().is_none());
    }
}
