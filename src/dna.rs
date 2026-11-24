//! DNA 序列工具（原 `gskits::dna`，随 gskits 依赖移除一并移植）。

/// 互补碱基表。未知的字节映射为 0（与原实现一致）。
static COMPLEMENT_TABLE: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'C' as usize] = b'G';
    table[b'G' as usize] = b'C';
    table[b'-' as usize] = b'-';
    table[b'*' as usize] = b'*';
    table[b'N' as usize] = b'N';
    table[b'a' as usize] = b't';
    table[b't' as usize] = b'a';
    table[b'c' as usize] = b'g';
    table[b'g' as usize] = b'c';
    table[b'n' as usize] = b'n';
    table
};

/// 返回反向互补序列（bytes in / bytes out）。
pub fn reverse_complement(dna: &[u8]) -> Vec<u8> {
    dna.iter().rev().map(|&b| COMPLEMENT_TABLE[b as usize]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ATCGTAGC"), b"GCTACGAT");
        assert_eq!(reverse_complement(b"-*"), b"*-");
        assert_eq!(reverse_complement(b"Nn"), b"nN");
    }
}
