# Shift-Or (Bitap) Algorithm - Technical Explanation

## Table of Contents
1. [Introduction](#introduction)
2. [Exact Matching Algorithm](#exact-matching-algorithm)
3. [Approximate Matching Algorithm](#approximate-matching-algorithm)
4. [Extended Multi-Word Implementation](#extended-multi-word-implementation)
5. [Complexity Analysis](#complexity-analysis)
6. [Implementation Details](#implementation-details)

## Introduction

The Shift-Or algorithm (also called Bitap or Baeza-Yates-Gonnet algorithm) is a **bit-parallel string matching** algorithm that uses bitwise operations to achieve fast pattern matching. It was originally invented by Bálint Dömölki in 1964 for exact matching and later extended by Ricardo Baeza-Yates and Gaston Gonnet in 1989.

### Key Innovation

Instead of comparing characters one by one, Shift-Or:
- Encodes the pattern as **bitmasks**
- Maintains a **state vector** using bitwise operations
- Processes text in **O(n)** time with **constant space**

### Why It's Fast for DNA

DNA sequences have a small alphabet (A, C, G, T, N), making bitmask operations extremely efficient:
- Only 5 bitmasks needed
- Each character lookup is O(1)
- Bit operations are hardware-accelerated

---

## Exact Matching Algorithm

### Core Idea

For a pattern P of length m:
1. Build bitmask B[c] for each character c
2. Maintain state vector D (initially all 1s)
3. For each text character, update D using bit operations
4. Match occurs when bit (m-1) becomes 0

### Bitmask Construction

For each character c in the alphabet:
```
B[c] = 111...111  (m bits, all set to 1)

For i = 0 to m-1:
    if P[i] == c:
        B[c][i] = 0  (clear bit i)
```

**Example**: Pattern "ACGT"
```
B['A'] = 1110  (bit 0 is 0, others are 1)
B['C'] = 1101  (bit 1 is 0)
B['G'] = 1011  (bit 2 is 0)
B['T'] = 0111  (bit 3 is 0)
```

### Search Algorithm

```
D = 111...111  (all bits set to 1)
match_bit = 1 << (m-1)  (bit at position m-1)

For each character text[j]:
    1. Get B[text[j]]
    2. D = ((D << 1) | 1) & B[text[j]]
    3. if (D & match_bit) == 0:
          Match found at position (j - m + 1)
```

### Why This Works

- **D[i] = 0** means "pattern[0...i] matches text ending at current position"
- Left shift prepares for next character
- OR with 1 sets bit 0 (pattern always starts matching)
- AND with B[] checks if current character matches pattern

**Time Complexity**: O(n) - one pass through text  
**Space Complexity**: O(|Σ|) - store bitmasks for alphabet

---

## Approximate Matching Algorithm

### Core Idea

Allow up to k errors (substitutions, insertions, deletions) by maintaining **k+1 state vectors**:
- D₀: exact matches (0 errors)
- D₁: matches with ≤1 error
- D₂: matches with ≤2 errors
- ...
- Dₖ: matches with ≤k errors

### Update Rules

For each text character:

**D₀ (exact matching)**:
```
D₀ = ((D₀ << 1) | 1) & B[text[j]]
```

**Dᵢ for i = 1 to k**:
```
substitution = D[i-1] << 1          # Previous error level, advance
insertion = D[i]                    # Current error level, don't advance
deletion = D[i-1]                   # Previous error level, don't advance
match = D[i] << 1                   # Current error level, advance

D[i] = ((substitution | insertion | deletion | match) | 1) & B[text[j]]
```

### Operations Explained

- **Substitution**: Accept mismatch, advance pattern (from Dᵢ₋₁ → Dᵢ)
- **Insertion**: Character in text not in pattern, advance text only
- **Deletion**: Character in pattern not in text, advance pattern only
- **Match**: Exact character match, advance both

### Detecting Matches

Check each Dᵢ from i=0 to k:
```
if (D[i] & match_bit) == 0:
    Match with at most i errors at position (j - m + 1)
    Report (position, i) and break  # Report lowest error level
```

**Time Complexity**: O(k × n) - k state vectors updated per character  
**Space Complexity**: O(k × |Σ|) - store k+1 state vectors

---

## Extended Multi-Word Implementation

### Problem

Standard Shift-Or uses a single 64-bit integer. For patterns > 64 bp, we need **multiple words**.

### Solution

Split pattern into ⌈m/64⌉ words:
- Word 0: pattern[0:64]
- Word 1: pattern[64:128]
- Word 2: pattern[128:192]
- ...

### Multi-Word Bitmasks

For each character c:
```
B[c] = [B[c][0], B[c][1], ..., B[c][w-1]]
```

where w = ⌈m/64⌉ words

### Multi-Word State Update

Maintain w state words: D[0], D[1], ..., D[w-1]

For each text character:
```
carry = 0
for word_idx = 0 to w-1:
    shifted = (D[word_idx] << 1) | carry

    # Extract carry for next word (bit 64)
    carry = 1 if (D[word_idx] & (1 << 63)) == 0 else 0

    # Apply bitmask
    D[word_idx] = (shifted | 1) & B[text[j]][word_idx]

# Check match in last word
if (D[w-1] & match_bit) == 0:
    Match found
```

### Carry Propagation

The **carry bit** propagates matching state across word boundaries:
- If bit 63 of word i is 0 (match), set carry=1 for word i+1
- This allows pattern matching to span multiple words

**Time Complexity**: O(w × n) where w = ⌈m/64⌉  
**Space Complexity**: O(w × |Σ|)

---

## Complexity Analysis

### Exact Matching (≤64 bp)

| Operation | Time | Space |
|-----------|------|-------|
| Preprocessing | O(m × |Σ|) | O(|Σ|) |
| Search | O(n) | O(1) |
| **Total** | **O(m×|Σ| + n)** | **O(|Σ|)** |

For DNA: |Σ| = 5, so preprocessing is O(5m) ≈ O(m)

### Approximate Matching (k errors)

| Operation | Time | Space |
|-----------|------|-------|
| Preprocessing | O(m × |Σ|) | O(|Σ|) |
| Search | O(k × n) | O(k) |
| **Total** | **O(m×|Σ| + k×n)** | **O(k × |Σ|)** |

### Extended Matching (>64 bp)

| Operation | Time | Space |
|-----------|------|-------|
| Preprocessing | O(m × |Σ|) | O(⌈m/64⌉ × |Σ|) |
| Search | O(⌈m/64⌉ × n) | O(⌈m/64⌉) |
| **Total** | **O(m×|Σ| + (m/64)×n)** | **O((m/64) × |Σ|)** |

---

## Implementation Details

### Pattern Length Boundaries

Our implementation has three variants:

1. **Standard Exact/Approximate** (m ≤ 64):
   - Single 64-bit integer
   - Optimal performance
   - Hardware-accelerated bit operations

2. **Extended** (64 < m ≤ 800):
   - Multiple 64-bit integers
   - Slight performance overhead from carry propagation
   - Still faster than character-by-character algorithms

### Bit Operation Counting

Per character processed:

- **Exact**: ~3 operations (shift, OR, AND)
- **Approximate (k=1)**: ~11 operations (6 for D₀, 5 more for D₁)
- **Approximate (k=2)**: ~19 operations
- **Approximate (k=3)**: ~27 operations
- **Extended**: ~7w operations where w = number of words

### Memory Usage

For DNA sequences (|Σ| = 5):

- **Exact**: 5 × 8 bytes = 40 bytes (for bitmasks)
- **Approximate (k=3)**: 5 × 8 × 4 = 160 bytes (4 state vectors)
- **Extended (800 bp)**: 5 × 8 × 13 = 520 bytes (13 words)

Compare to KMP: O(m) bytes for LPS array

### Performance Characteristics

**Best Case**: Pattern occurs frequently  
**Worst Case**: Pattern never occurs (still O(n))  
**Average Case**: O(n) regardless of matches

The algorithm is **oblivious** - running time depends only on text length, not on matches found.

---

## Advantages

✅ **Linear time** O(n) for exact matching  
✅ **Simple to implement** - mostly bitwise operations  
✅ **Hardware-friendly** - CPUs optimize bit operations  
✅ **Predictable performance** - no worst-case slowdowns  
✅ **Naturally extends** to approximate matching  

## Disadvantages

❌ **Pattern length limit** - standard version limited to 64 bp  
❌ **Not cache-friendly** - random memory access for bitmasks  
❌ **Slower preprocessing** - O(m × |Σ|) vs O(m) for KMP  
❌ **Approximate matching overhead** - O(k × n) can be slower for large k  

---

## Comparison with Other Algorithms

| Algorithm | Preprocessing | Search | Space | Best For |
|-----------|--------------|--------|-------|----------|
| KMP | O(m) | O(n) | O(m) | Any pattern length |
| Boyer-Moore | O(m+|Σ|) | O(n/m) avg | O(m+|Σ|) | Long patterns, large alphabet |
| **Shift-Or** | **O(m×|Σ|)** | **O(n)** | **O(|Σ|)** | **Short patterns, small alphabet** |

**Shift-Or excels for DNA** because:
- Small alphabet (|Σ| = 5)
- Patterns often ≤ 64 bp
- Bit operations are very fast
- Naturally handles approximate matching

---

## References

1. Dömölki, B. (1964). "A universal automaton for searching and identification." Studia Scientiarum Mathematicarum Hungarica.

2. Baeza-Yates, R., & Gonnet, G. H. (1992). "A new approach to text searching." Communications of the ACM, 35(10), 74-82.

3. Wu, S., & Manber, U. (1992). "Fast text searching: allowing errors." Communications of the ACM, 35(10), 83-91.

4. Navarro, G., & Raffinot, M. (2002). "Flexible Pattern Matching in Strings." Cambridge University Press.
