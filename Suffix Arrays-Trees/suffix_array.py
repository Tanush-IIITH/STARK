"""Manber-Myers suffix array construction.

This module provides a single helper `manber_myers_suffix_array` that
materializes the suffix array for an arbitrary string.  The implementation
closely follows the classic Manber-Myers doubling strategy, carrying along
explicit commentary so that each stage of the algorithm remains easy to trace.
"""

from __future__ import annotations

import time
from typing import Callable, Iterable, List


def brute_force_suffix_array(text: str) -> List[int]:
	"""Return the suffix array of *text* via the naive O(n^2 log n) routine.

	This mirrors the conceptual steps one would take manually: enumerate all
	suffixes, sort them lexicographically, and report the starting offsets of
	the sorted suffixes.  Although prohibitively slow for large inputs, the
	routine is helpful for educational purposes and for validating more
	sophisticated implementations.
	"""

	n = len(text)
	if n == 0:
		return []

	# Materialize every suffix together with its origin index.
	suffixes = [(text[i:], i) for i in range(n)]

	# Sorting by the suffix value itself matches the lexicographic ordering.
	suffixes.sort(key=lambda item: item[0])

	# Strip away the suffix strings, leaving only the original indices.
	return [index for _, index in suffixes]


def manber_myers_suffix_array(text: str) -> List[int]:
	"""Return the suffix array of *text* using the Manber-Myers algorithm.

	The algorithm iteratively doubles the length of the compared prefix.  At
	each stage we maintain a `rank` array that stores the equivalence class of
	every suffix with respect to the first ``2^k`` characters that were already
	processed.  Sorting suffix indices by the pair of ranks ``(rank[i],
	rank[i + k])`` yields the correct order for ``2^(k+1)`` characters; the
	new order allows us to compute the next set of ranks and continue doubling
	until all suffixes are distinguishable.

	Args:
		text: Input string for which the suffix array should be produced.

	Returns:
		A list whose *i*-th element records the starting index of the *i*-th
		lexicographically smallest suffix of *text*.
	"""

	n = len(text)
	if n == 0:
		return []

	# Initial ordering and ranks rely on single characters.
	suffix_array = list(range(n))
	rank = [ord(ch) for ch in text]

	k = 1
	def counting_sort(indices: Iterable[int], key_fn: Callable[[int], int], key_range: int) -> List[int]:
		"""Stable counting sort specialised for suffix indices."""

		indices_list = list(indices)
		counts = [0] * key_range
		for idx in indices_list:
			counts[key_fn(idx)] += 1

		for i in range(1, key_range):
			counts[i] += counts[i - 1]

		output = [0] * len(indices_list)
		for idx in reversed(indices_list):
			key = key_fn(idx)
			counts[key] -= 1
			output[counts[key]] = idx
		return output

	while k < n:
		max_rank = max(rank)

		# Sort by the second half (rank at idx + k).  The sentinel -1 is mapped to 0.
		second_key_range = max_rank + 2
		second_key = lambda idx: (rank[idx + k] + 1) if idx + k < n else 0
		suffix_array = counting_sort(suffix_array, second_key, second_key_range)

		# Sort by the first half (current rank).
		first_key_range = max_rank + 1
		first_key = lambda idx: rank[idx]
		suffix_array = counting_sort(suffix_array, first_key, first_key_range)

		new_rank = [0] * n
		new_rank[suffix_array[0]] = 0
		distinct_ranks = 0

		# Walk the sorted suffix indices, emitting a new rank whenever the
		# neighbouring (rank, next_rank) pair changes.
		for pos in range(1, n):
			prev_idx = suffix_array[pos - 1]
			curr_idx = suffix_array[pos]

			prev_pair = (rank[prev_idx], rank[prev_idx + k] if prev_idx + k < n else -1)
			curr_pair = (rank[curr_idx], rank[curr_idx + k] if curr_idx + k < n else -1)

			if curr_pair != prev_pair:
				distinct_ranks += 1

			new_rank[curr_idx] = distinct_ranks

		rank = new_rank

		# Once all ranks are unique we can short-circuit; further doubling
		# would not change the ordering.
		if distinct_ranks == n - 1:
			break

		k <<= 1

	return suffix_array


def locate_pattern(text: str, pattern: str, suffix_array: List[int]) -> List[int]:
	"""Return sorted starting indices where *pattern* occurs in *text*."""

	if not pattern:
		return list(range(len(text)))

	def compare(index: int) -> int:
		probe = text[index : index + len(pattern)]
		if probe == pattern:
			return 0
		return -1 if probe < pattern else 1

	# Locate the first suffix whose prefix is >= pattern.
	left, right = 0, len(suffix_array)
	while left < right:
		mid = (left + right) // 2
		cmp = compare(suffix_array[mid])
		if cmp < 0:
			left = mid + 1
		else:
			right = mid
	first = left

	# Locate the first suffix whose prefix is > pattern.
	left, right = first, len(suffix_array)
	while left < right:
		mid = (left + right) // 2
		cmp = compare(suffix_array[mid])
		if cmp > 0:
			right = mid
		else:
			left = mid + 1
	last = left

	return [suffix_array[i] for i in range(first, last) if text[suffix_array[i] :].startswith(pattern)]


def _interactive_cli() -> None:
	"""Drive a simple REPL for choosing algorithms and querying patterns."""

	algorithms: dict[str, Callable[[str], List[int]]] = {
		"brute": brute_force_suffix_array,
		"manber": manber_myers_suffix_array,
		"manber-myers": manber_myers_suffix_array,
		"naive": brute_force_suffix_array,
	}

	print("Available algorithms: brute, manber")
	while True:
		choice = input("Select algorithm: ").strip().lower()
		if choice in algorithms:
			break
		print("Unrecognised choice. Please enter 'brute' or 'manber'.")

	text = input("Enter the text: ")
	pattern = input("Enter the pattern to search: ")
	show_choice = input("Show suffix array details? [y/N]: ").strip().lower()
	show_details = show_choice in {"y", "yes"}

	builder = algorithms[choice]
	start_time = time.perf_counter()
	suffix_array = builder(text)
	elapsed = time.perf_counter() - start_time

	print()
	print(f"Suffix array using '{choice}' ({elapsed:.6f} seconds):")
	if show_details:
		if not suffix_array:
			print("  [empty text]")
		else:
			print("  order | start | suffix")
			for order, idx in enumerate(suffix_array):
				suffix = text[idx:]
				print(f"  {order:5d} | {idx:5d} | {suffix}")
	else:
		print("  (suffix array output suppressed by user request)")

	occurrences = locate_pattern(text, pattern, suffix_array)
	print()
	print(f"Pattern occurrences ({len(occurrences)} matches): {occurrences}")


if __name__ == "__main__":
	try:
		_interactive_cli()
	except KeyboardInterrupt:
		print("\nInterrupted by user.")


__all__ = [
	"brute_force_suffix_array",
	"manber_myers_suffix_array",
	"locate_pattern",
]
