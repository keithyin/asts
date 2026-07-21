use std::{fmt::Debug, str::FromStr};

use serde::Serialize;

/// parse the range str to Vector of tuple
/// range fmt: b1:e1,b2:e2,be3,b4:e4. corresponds to [[b1, e1], [b2, e2], [be3, be3], [b4, e4]]. left inclusive, right inclusive
#[allow(unused)]
pub fn range_parser<T>(range_str: &str) -> Vec<(T, T)>
where
    T: FromStr,
    <T as FromStr>::Err: Debug,
{
    range_str
        .trim()
        .split(",")
        .map(|single_range| {
            let items = single_range
                .trim()
                .split(":")
                .filter(|item| item.trim().len() != 0)
                .collect::<Vec<_>>();
            return if items.len() == 2 {
                (
                    items[0].parse::<T>().expect(&format!("{}", items[0])),
                    items[1].parse::<T>().expect(&format!("{}", items[1])),
                )
            } else {
                (
                    items[0].parse::<T>().expect(&format!("{}", items[0])),
                    items[0].parse::<T>().expect(&format!("{}", items[0])),
                )
            };
        })
        .collect::<Vec<_>>()
}

#[derive(Debug, Clone, Serialize)]
pub struct Range<T> {
    range: Vec<(T, T)>,
}

impl<T> Range<T>
where
    T: PartialOrd,
    T: FromStr,
    <T as FromStr>::Err: Debug,
{
    pub fn new(range_str: &str) -> Self {
        Self {
            range: range_parser::<T>(range_str),
        }
    }

    pub fn new_range(range: Vec<(T, T)>) -> Self {
        Self { range }
    }

    pub fn within_range(&self, v: T) -> bool {
        let mut result = false;
        for (b, e) in &self.range {
            if *b <= v && v <= *e {
                result = true;
                break;
            }
        }
        result
    }
}

impl Range<usize> {
    pub fn shift(&mut self, shift: i64) {
        self.range.iter_mut().for_each(|(start, end)| {
            let new_start = (*start as i64 + shift).max(0) as usize;
            let new_end = (*end as i64 + shift).max(0) as usize;
            *start = new_start;
            *end = new_end;
        });
    }
}

impl Range<i64> {
    pub fn shift(&mut self, shift: i64) {
        self.range.iter_mut().for_each(|(start, end)| {
            let new_start = *start as i64 + shift;
            let new_end = *end as i64 + shift;
            *start = new_start;
            *end = new_end;
        });
    }
}
