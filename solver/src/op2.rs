//! Minimal Nastran OP2: one `OUGV1` displacement table (SORT1, real).
//!
//! Layout follows pyNastran's writer (little-endian Fortran records): a
//! NASTRAN tape header, the `OUGV1` table header, a 584-byte table 3 and an
//! 8-wide grid record (`nid*10+device`, grid type, T1..R3).

use crate::analysis::SolveOutput;
use crate::model::Model;

const DEVICE: i32 = 1;

pub fn write_op2(model: &Model, out: &SolveOutput) -> Vec<u8> {
    let mut buf = Vec::new();
    write_tape_header(&mut buf);
    let blocks = displacement_blocks(model, out);
    if blocks.is_empty() {
        push_i32s(&mut buf, &[4, 0, 4]);
        return buf;
    }
    write_table_header(&mut buf, "OUGV1   ", "OUG1    ");
    let mut itable = -3;
    for (i, block) in blocks.iter().enumerate() {
        write_table3(&mut buf, block, itable, i == 0);
        itable -= 1;
        write_table4(&mut buf, model, block, itable);
        itable -= 1;
    }
    push_i32s(&mut buf, &[4, itable, 4, 4, 1, 4, 4, 0, 4]);
    push_i32s(&mut buf, &[4, 0, 4]);
    push_i32s(&mut buf, &[4, 0, 4]);
    buf
}

struct Block {
    analysis: i32,
    isubcase: i32,
    field5_i: i32,
    field5_f: Option<f32>,
    field6: f32,
    field7: f32,
    u: Vec<[f64; 3]>,
    ur: Vec<[f64; 3]>,
}

fn displacement_blocks(model: &Model, out: &SolveOutput) -> Vec<Block> {
    let mut all = Vec::new();
    if out.cases.is_empty() {
        all.extend(expand(
            1,
            &out.u,
            &out.ur,
            &out.frequencies,
            &out.buckles,
            &out.modes,
            &out.mode_ur,
        ));
    } else {
        for (i, case) in out.cases.iter().enumerate() {
            all.extend(expand(
                i as i32 + 1,
                &case.u,
                &case.ur,
                &case.frequencies,
                &case.buckles,
                &case.modes,
                &case.mode_ur,
            ));
        }
    }
    if all.is_empty() {
        all.push(static_block(1, &out.u, &out.ur));
    }
    let _ = model;
    all
}

fn expand(
    isub: i32,
    u: &[[f64; 3]],
    ur: &[[f64; 3]],
    freq: &[f64],
    buck: &[f64],
    modes: &[Vec<[f64; 3]>],
    mode_ur: &[Vec<[f64; 3]>],
) -> Vec<Block> {
    if !modes.is_empty() && !buck.is_empty() {
        return modes
            .iter()
            .enumerate()
            .map(|(i, uu)| {
                buckle_block(
                    isub,
                    (i + 1) as i32,
                    buck.get(i).copied().unwrap_or(0.0),
                    uu,
                    mode_ur.get(i).map(|v| v.as_slice()).unwrap_or(&[]),
                )
            })
            .collect();
    }
    if !modes.is_empty() && (!freq.is_empty() || modes.len() > 1) {
        return modes
            .iter()
            .enumerate()
            .map(|(i, uu)| {
                let f = freq.get(i).copied().unwrap_or(0.0);
                let w2 = if f > 0.0 {
                    (2.0 * std::f64::consts::PI * f).powi(2) as f32
                } else {
                    0.0
                };
                Block {
                    analysis: 2,
                    isubcase: isub,
                    field5_i: (i + 1) as i32,
                    field5_f: None,
                    field6: w2,
                    field7: f as f32,
                    u: uu.clone(),
                    ur: mode_ur.get(i).cloned().unwrap_or_default(),
                }
            })
            .collect();
    }
    vec![static_block(isub, u, ur)]
}

fn static_block(isub: i32, u: &[[f64; 3]], ur: &[[f64; 3]]) -> Block {
    Block {
        analysis: 1,
        isubcase: isub,
        field5_i: 0,
        field5_f: None,
        field6: 0.0,
        field7: 0.0,
        u: u.to_vec(),
        ur: ur.to_vec(),
    }
}

fn buckle_block(isub: i32, mode: i32, factor: f64, u: &[[f64; 3]], ur: &[[f64; 3]]) -> Block {
    Block {
        analysis: 8,
        isubcase: isub,
        field5_i: mode,
        field5_f: None,
        field6: factor as f32,
        field7: 0.0,
        u: u.to_vec(),
        ur: ur.to_vec(),
    }
}

fn write_tape_header(buf: &mut Vec<u8>) {
    push_i32s(buf, &[4, 3, 4]);
    let tape = b"NASTRAN FORT TAPE ID CODE - ";
    // MSC date record: 12, day, month, year-2000, then markers.
    let mut rec = Vec::new();
    push_i32s(&mut rec, &[12, 25, 9, 26, 12, 4, 7, 4, 28]);
    rec.extend_from_slice(tape);
    push_i32s(&mut rec, &[28]);
    buf.extend_from_slice(&rec);
    push_i32s(buf, &[4, 2, 4, 8]);
    buf.extend_from_slice(b"XXXXXXXX");
    push_i32s(buf, &[8, 4, -1, 4, 4, 0, 4]);
}

fn write_table_header(buf: &mut Vec<u8>, table: &str, sub: &str) {
    push_i32s(buf, &[4, 2, 4, 8]);
    let mut name = [b' '; 8];
    let tb = table.as_bytes();
    name[..tb.len().min(8)].copy_from_slice(&tb[..tb.len().min(8)]);
    buf.extend_from_slice(&name);
    push_i32s(buf, &[8, 4, -1, 4, 4, 7, 4, 28, 102, 0, 0, 0, 512, 0, 0, 28]);
    push_i32s(buf, &[4, -2, 4, 4, 1, 4, 4, 0, 4, 4, 7, 4, 28]);
    let mut subn = [b' '; 8];
    let sb = sub.as_bytes();
    subn[..sb.len().min(8)].copy_from_slice(&sb[..sb.len().min(8)]);
    buf.extend_from_slice(&subn);
    // month, day, year-2000, 0, 1, trailing record length
    push_i32s(buf, &[9, 25, 26, 0, 1, 28]);
}

fn write_table3(buf: &mut Vec<u8>, block: &Block, itable: i32, first: bool) {
    if first {
        push_i32s(buf, &[4, itable, 4, 4, 1, 4, 4, 0, 4, 4, 146, 4]);
    } else {
        push_i32s(buf, &[4, 146, 4]);
    }
    let approach = block.analysis * 10 + DEVICE;
    let mut words = [0i32; 50];
    words[0] = approach;
    words[1] = 1; // OUG displacement
    words[3] = block.isubcase;
    if let Some(f) = block.field5_f {
        words[4] = f.to_bits() as i32;
    } else {
        words[4] = block.field5_i;
    }
    if block.analysis == 2 || block.analysis == 8 {
        words[5] = block.field6.to_bits() as i32;
    }
    if block.analysis == 2 {
        words[6] = block.field7.to_bits() as i32;
    }
    words[8] = 1; // real
    words[9] = 8; // num_wide
    // words[22] thermal = 0
    push_i32(buf, 584);
    for w in words {
        push_i32(buf, w);
    }
    for _ in 0..3 {
        buf.extend_from_slice(&[b' '; 128]);
    }
    push_i32(buf, 584);
}

fn write_table4(buf: &mut Vec<u8>, model: &Model, block: &Block, itable: i32) {
    let n = model.node_ids.len().min(block.u.len());
    let ntotal = (n as i32) * 8;
    let nbytes = ntotal * 4;
    push_i32s(buf, &[4, itable, 4, 4, 1, 4, 4, 0, 4, 4, ntotal, 4, nbytes]);
    for i in 0..n {
        let nid = model.node_ids[i];
        push_i32(buf, nid.saturating_mul(10) + DEVICE);
        push_i32(buf, 1); // GRID
        let u = block.u[i];
        let r = block.ur.get(i).copied().unwrap_or([0.0; 3]);
        for v in [u[0], u[1], u[2], r[0], r[1], r[2]] {
            push_f32(buf, v as f32);
        }
    }
    push_i32(buf, nbytes);
}

fn push_i32(buf: &mut Vec<u8>, v: i32) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn push_f32(buf: &mut Vec<u8>, v: f32) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn push_i32s(buf: &mut Vec<u8>, vs: &[i32]) {
    for v in vs {
        push_i32(buf, *v);
    }
}

/// Read the first OUGV1 grid record back. Returns (nid, t1, t2, t3).
pub fn read_ougv1_first(bytes: &[u8]) -> Option<(i32, [f32; 6])> {
    let table = b"OUGV1";
    let pos = bytes.windows(8).position(|w| &w[..5] == table)?;
    // Data follows the 584-byte table-3 record. Find 584 marker after the name.
    let mut i = pos;
    let mut saw_584 = false;
    while i + 8 <= bytes.len() {
        let n = i32::from_le_bytes(bytes[i..i + 4].try_into().ok()?);
        if n == 584 && !saw_584 {
            saw_584 = true;
            i += 4 + 584;
            if i + 4 <= bytes.len() {
                let end = i32::from_le_bytes(bytes[i..i + 4].try_into().ok()?);
                if end == 584 {
                    i += 4;
                }
            }
            continue;
        }
        if saw_584 && n > 8 && n % 4 == 0 && i + 4 + 32 <= bytes.len() {
            // Could be the ntotal marker (number of words) not the byte length.
            let words = n;
            if words % 8 == 0 && words < 1_000_000 {
                // next pattern in our writer: [4, ntotal, 4, nbytes, data...]
                // We landed on ntotal only if the previous ints were markers.
            }
        }
        i += 4;
        if saw_584 && n > 32 && n % 32 == 0 && i + n as usize <= bytes.len() {
            let data = &bytes[i..i + n as usize];
            let nid_dev = i32::from_le_bytes(data[0..4].try_into().ok()?);
            let mut v = [0.0; 6];
            for k in 0..6 {
                let o = 8 + k * 4;
                v[k] = f32::from_le_bytes(data[o..o + 4].try_into().ok()?);
            }
            return Some((nid_dev / 10, v));
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ougv1_contains_table_name() {
        let model = crate::parse_model(
            "SOL 101\nCEND\nSPC=1\nLOAD=1\nBEGIN BULK\n\
             GRID,1,,0.,0.,0.\nGRID,2,,1.,0.,0.\n\
             CROD,1,1,1,2\nPROD,1,1,2.\nMAT1,1,100.,,0.\n\
             FORCE,1,2,,10.,1.,0.,0.\nSPC1,1,123456,1\nENDDATA\n",
        )
        .unwrap();
        let out = crate::analysis::solve(model.clone()).unwrap();
        let bytes = write_op2(&model, &out);
        assert!(bytes.windows(8).any(|w| &w[..5] == b"OUGV1"), "no OUGV1");
        let (nid, v) = read_ougv1_first(&bytes).expect("grid record");
        assert_eq!(nid, model.node_ids[0]);
        assert!((v[0] as f64 - out.u[0][0]).abs() < 1e-5, "t1 {} vs {}", v[0], out.u[0][0]);
        let tip = 0.05f32.to_le_bytes();
        assert!(bytes.windows(4).any(|w| w == tip), "OP2 missing tip ux 0.05");
    }
}
