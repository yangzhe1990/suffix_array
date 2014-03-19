#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>

/* temporary don't implement */
#define _suffix_array_malloc malloc

#define MASK_N 0x40000000L
// MAXN < MASK_N. MASK_N must be a power of 2 and MASK_N - 1 +
// MASK_N must not overflow a signed int

struct timeval tv;

/* in the _suffix_array functions, padding $ is not included in either
 * the pattern or the alphabet */
/* type = 0, 1, stand for S, L, respectively. The type of the last char is L */
/* temporary don't touch lcp */

template <typename T>
int _suffix_array_sort_substring_s(T *pattern, int length, int new_length,
				   int alphabet_size,
				   int *orig_pos, int *sorted_suffix,
				   int *bucket_end, int *bucket_p,
				   int *sorted_substrings,
				   int **p_new_pattern,
				   int **p_new_bucket_end,
				   int *p_new_r_alphabet_size,
				   int *p_new_r_length) {
	/* fill the type S positions into the bucket */
	memcpy(bucket_p, bucket_end, alphabet_size * sizeof(int));
	//memset(sorted_suffix, 0, length * sizeof(int));
	int i;
	int j = 0;
	int same = 0; // same means the same char
	T *c = pattern;
	int c1 = c[0], c0 = c1, this_chr; // we should avoid overflow of signed char
	int *b = sorted_suffix + bucket_p[c0];
	int r = 0;
	// assume -1 is a type S position. That means, if 0 is a type
	// S position, it won't enter the forward update branch.
	for (i = 1; i < length; ++i) {
		this_chr = *(++c);
		if (c1 < this_chr) {
			if (c1 != c0) {
				bucket_p[c0] = b - sorted_suffix;
				b = sorted_suffix + bucket_p[c0 = c1];
			}

			--b;
			if (r == 1) {
				// same - 1 is a type L position, mark same.
				orig_pos[j] = ~same;
				*b = ~same ^ MASK_N;
				r = 0;
			}
			else {
				orig_pos[j] = same;
				*b = same ^ MASK_N;
			}
			++j;

			for (++same; same < i; ++same) {
				orig_pos[j++] = same;
				*(--b) = same ^ MASK_N;
			}
		}
		else if (c1 > this_chr) {
			same = i;
			r = 1;
		}
		c1 = this_chr;
	}
	// if same != length - 1, then the last same block is consist
	// of all type L positions.
	bucket_p[c0] = b - sorted_suffix;
	// alphabet_size is at least 2
	for (i = alphabet_size - 1; i >= 0; --i)
		if ((j = bucket_p[i]) < bucket_end[i])
			sorted_suffix[j] ^= MASK_N;
	/* no need to assign equal value for type L suffixes and no
	 * need to set the order for type L suffixes. */

	/* sort the substrings */
	int *s_bucket_begin =
		(int *)_suffix_array_malloc(alphabet_size * sizeof(int));
	// In fact we can reuse bucket_p, but how to compute the
	// s_bucket_begin in O(1) space? by +- sign of bucket_p?

	j = 0;
	for (i = 0; i < alphabet_size; ++i) {
		s_bucket_begin[i] = j;
		j += bucket_end[i] - bucket_p[i];
	}

	int *last_deleg =
		(int *)_suffix_array_malloc(alphabet_size * sizeof(int));
	memset(last_deleg, -1, alphabet_size * sizeof(int));

	// init bucket head for filling type L positions
	bucket_p[0] = 0;
	memcpy(bucket_p + 1, bucket_end, (alphabet_size - 1) * sizeof(int));

	assert(length > 1);
	// pattern[length] is $, which is the smallest char. The
	// position length - 1 is a type L position.
	i = length - 1;
	c0 = pattern[i];
	b = sorted_suffix + bucket_p[c0];
	// if i - 1 is type L pos, mark i;
	if (pattern[i - 1] >= c0)
		*b = ~i;
	else
		*b = i;
	++b;
	int s_c0 = 0;
	int *s_b = sorted_substrings + s_bucket_begin[s_c0];
	int s;
	int equal_deleg = 0;
	int new_alphabet_size = new_length;
	for (s = 0; s < length; ++s) {
		i = sorted_suffix[s];
		if (i < 0) {
			i ^= -1;
			if (i >= MASK_N)
				i ^= MASK_N;
			else
				equal_deleg = s;
			--i;
			// update forward
			if ((c1 = pattern[i]) != c0) {
				bucket_p[c0] = b - sorted_suffix;
				b = sorted_suffix + bucket_p[c0 = c1];
			}
			// if i - 1 is type L pos, mark i;
			if (i > 0 && pattern[i - 1] >= c0)
				i ^= -1;

			if (equal_deleg == last_deleg[c0])
				*b = i ^ MASK_N;
			else {
				*b = i;
				last_deleg[c0] = equal_deleg;
			}
			++b;
		}
		else {
			if (i >= MASK_N)
				i ^= MASK_N;
			else
				equal_deleg = s;
			if (i == 0)
				continue;
			// update backward to determine a substring's rank
			--i;
			if ((c1 = pattern[i]) != s_c0) {
				s_bucket_begin[s_c0] = s_b - sorted_substrings;
				s_b = sorted_substrings
					+ s_bucket_begin[s_c0 = c1];
			}

			// if i - 1 is type L pos, mark i;
			if (i > 0 && pattern[i - 1] > s_c0)
				i ^= -1;

			if (equal_deleg == last_deleg[s_c0]) {
				*s_b = i ^ MASK_N;
				--new_alphabet_size;
			}
			else {
				*s_b = i;
				last_deleg[s_c0] = equal_deleg;
			}
			++s_b;
		}
	}
	free(s_bucket_begin);
	free(last_deleg);

	/* compute the new pattern for recursion */
	gettimeofday(&tv, NULL);
	fprintf(stderr, " %d prepare rec at %d.%06d\n",
		length, tv.tv_sec, tv.tv_usec);

	if (new_alphabet_size > 1 && new_alphabet_size < new_length) {
		// save the rank of substrings of each pos in
		// sorted_suffix[pos]. If a position has no tie, its
		// rank ^ -1 will be saved.
		int k = sorted_substrings[0];
		if (k < 0)
			k ^= -1;
		//assert(k < MASK_N);
		j = 0;
		for (s = 1; s < new_length;) {
			i = k;

			do {
				k = sorted_substrings[s];
				if (k < 0)
					k ^= -1;
				if (k >= MASK_N)
					sorted_suffix[k ^ MASK_N] = j;
				else
					break;
				++s;
			} while (s < new_length);
			if (j + 1 == s)
				sorted_suffix[i] = ~j;
			else
				sorted_suffix[i] = j;

			j = s++;
		}
		if (j < new_length) // == new_length - 1
			sorted_suffix[k] = ~j;

		// reduce the new pattern
		int new_r_length = 0;
		for (i = 0; i < new_length; ) {
			s = orig_pos[i];
			if (s < 0)
				k = sorted_suffix[~s];
			else
				k = sorted_suffix[s];
			if (k >= 0) {
				if (s < 0)
					sorted_suffix[~s] = new_r_length;
				else
					sorted_suffix[s] = new_r_length;
				orig_pos[new_r_length++] = s;
				// any char following a char appeared
				// many times in the new pattern
				// should be kept.
				for (++i; (k >= 0) && i < new_length; ++i) {
					s = orig_pos[i];
					if (s < 0) {
						k = sorted_suffix[~s];
						sorted_suffix[~s] = new_r_length;
						// save the position
						// in the new pattern
					}
					else {
						k = sorted_suffix[s];
						sorted_suffix[s] = new_r_length;
					}
					orig_pos[new_r_length++] = s;
				}
			}
			else {
				// mark sorted_suffix[s] to notice that s is removed
				if (s < 0)
					sorted_suffix[~s] = -1;
				else
					sorted_suffix[s] = -1;
				++i;
			}
		}

		// generate the new alphabet and the new pattern.
		int *new_pattern =
			(int *)_suffix_array_malloc(
				new_length * sizeof(int));
			// sorted_suffix; // we can reuse it as new_pattern
			// FIXME. remove test code;
		int *new_bucket_end =
			/*
			(int *)_suffix_array_malloc(
				new_alphabet_size * sizeof(int));
			*/
			sorted_substrings; // we can reuse it as new_bucket_end;

		int r_i = new_r_length;
		// find the first pos which is in the new pattern. it
		// exists.
		j = 0;
		for (s = 0; s < new_length; ++s) {
			k = sorted_substrings[s];
			if (k < 0) {
				//assert(k > ~MASK_N);
				i = sorted_suffix[~k];
			}
			else {
				//assert(k < MASK_N);
				i = sorted_suffix[k];
			}
			if (i >= 0)
				break;
			else {
				// this pos is not in the new pattern
				new_pattern[r_i] = j;
				orig_pos[r_i] = k;
				++r_i;
			}
		}
		// this pos is in the new pattern
		new_pattern[i] = j;
		new_r_length = 1;

		for (++s; s < new_length; ++s) {
			// the first pos after a in-the-new-pattern pos
			k = sorted_substrings[s];
			if (k < 0)
				if (k <= ~MASK_N) {
					i = sorted_suffix[~k ^ MASK_N];
					//assert(i >= 0);
					goto end;
				}
				else
					i = sorted_suffix[~k];
			else
				if (k >= MASK_N) {
					i = sorted_suffix[k ^ MASK_N];
					//assert(i >= 0);
					goto end;
				}
				else
					i = sorted_suffix[k];
			new_bucket_end[j++] = new_r_length;
			if (i < 0) {
				// this pos is not in the new pattern
				new_pattern[r_i] = j;
				orig_pos[r_i] = k;
				++r_i;
				for (++s; s < new_length; ++s) {
					k = sorted_substrings[s];
					if (k < 0) {
						//assert(k > ~MASK_N);
						i = sorted_suffix[~k];
					}
					else {
						//assert(k < MASK_N);
						i = sorted_suffix[k];
					}
					if (i >= 0)
						break;
					else {
						// this pos is not in
						// the new pattern
						new_pattern[r_i] = j;
						orig_pos[r_i] = k;
						++r_i;
					}
				}
				if (s == new_length)
					break;
			}
		end:
			// the current pos is in the new pattern
			new_pattern[i] = j;
			++new_r_length;
		}
		new_bucket_end[j++] = new_r_length;

		*p_new_r_alphabet_size = j;
		*p_new_r_length = new_r_length;
		*p_new_pattern = new_pattern;
		*p_new_bucket_end = new_bucket_end;
		return new_alphabet_size;
	}
	else
		return 0;
}

template <typename T>
void _suffix_array_post_recursion_s(T *pattern, int length,
				    int alphabet_size,
				    int *sorted_suffix,
				    int *bucket_end, int *bucket_p) {
	// init bucket head for filling type L positions
	bucket_p[0] = 0;
	memcpy(bucket_p + 1, bucket_end, (alphabet_size - 1) * sizeof(int));

	assert(length > 1);
	int i = length - 1;
	int c0 = pattern[i], c1;
	int *b = sorted_suffix + bucket_p[c0];
	// if i - 1 is a type L position, mark i
	if (pattern[i - 1] >= c0)
		*b = ~i;
	else
		*b = i;
	++b;
	int s;
	for (s = 0; s < length; ++s) {
		i = sorted_suffix[s];
		if (i < 0) {
			i ^= -1;
			sorted_suffix[s] = i;
			--i;
			if ((c1 = pattern[i]) != c0) {
				bucket_p[c0] = b - sorted_suffix;
				b = sorted_suffix + bucket_p[c0 = c1];
			}
			// if i - 1 is a type L position, mark i
			if (i > 0 && pattern[i - 1] >= c0)
				*b = ~i;
			else
				*b = i;
			++b;
		}
	}
}

void _suffix_array_init(unsigned char *pattern, int length, int *bucket_end) {
	gettimeofday(&tv, NULL);
	fprintf(stderr, "init %d at %d.%06d\n", length, tv.tv_sec, tv.tv_usec);

	int alphabet_size = 256;

	/* construct the bucket */
	// char_count must be inited to 0
	int *char_count = (int *)_suffix_array_malloc(alphabet_size * sizeof(int));
	memset(char_count, 0, alphabet_size * sizeof(int));
	unsigned char *p = pattern + length;
	unsigned char *c;
	for (c = pattern; c < p; ++c)
		++char_count[*c];
	int tmp_int1 = 0;
	int i;
	for (i = 0; i < alphabet_size; ++i) {
		/* tmp_int1 : bucket_end[i] */
		tmp_int1 += char_count[i];
		bucket_end[i] = tmp_int1;
	}
	free(char_count);
}

template <typename T>
int _suffix_array_sort_substring_l(T *pattern, int length, int new_length,
				   int alphabet_size,
				   int *orig_pos, int *sorted_suffix,
				   int *bucket_end, int *bucket_p,
				   int *sorted_substrings,
				   int **p_new_pattern,
				   int **p_new_bucket_end,
				   int *p_new_r_alphabet_size,
				   int *p_new_r_length) {
	/* fill the type L positions into the bucket */
	bucket_p[0] = 0;
	memcpy(bucket_p + 1, bucket_end, (alphabet_size - 1) * sizeof(int));
	//memset(sorted_suffix, 0, length * sizeof(int));
	int i;
	int j = new_length;
	T *c = pattern + length - 1;
	int c1 = c[0], c0 = c1, this_chr;
	int *b = sorted_suffix + bucket_p[c0];
	int r = 1;
	// assume -1 is a type L position. That means, if 0 is a type
	// L position, it won't enter the backward update branch.
	for (i = length - 1; i > 0; --i) {
		this_chr = *(--c);
		if (this_chr > c1) {
			// this_chr is L
			if (r == 1) {
				// the next chr is L, save it
				orig_pos[--j] = i;
				*(b++) = i ^ MASK_N;
				bucket_p[c0] = b - sorted_suffix;
				b = sorted_suffix + bucket_p[c0 = this_chr];
			}
			else {
				if (c1 != c0) {
					bucket_p[c0] = b - sorted_suffix;
					b = sorted_suffix
						+ bucket_p[c0 = this_chr];
				}
				r = 1;
			}
		}
		else if (this_chr < c1) {
			// this_chr is S
			if (r == 1) {
				// the next chr is L, save it, with mark
				orig_pos[--j] = ~i;
				*(b++) = ~i ^ MASK_N;
				r = 0;
			}
		}
		else {
			if (r == 1) {
				// the next chr is L, this chr is also
				// L, save the next chr
				orig_pos[--j] = i;
				*(b++) = i ^ MASK_N;
			}
		}
		c1 = this_chr;
	}
	// if the first chr is L, save it.
	if (r == 1) {
		orig_pos[--j] = i;
		*(b++) = i ^ MASK_N;
	}
	bucket_p[c0] = b - sorted_suffix;
	// alphabet_size is at least 1
	if ((j = bucket_p[0] - 1) >= 0)
		sorted_suffix[j] ^= MASK_N;
	for (i = 1; i < alphabet_size; ++i)
		if ((j = bucket_p[i] - 1) >= bucket_end[i - 1])
			sorted_suffix[j] ^= MASK_N;
	/* no need to assign equal value for type L suffixes and no
	 * need to set the order for type L suffixes. */

	/* sort the substrings */
	int *l_bucket_end =
		(int *)_suffix_array_malloc(alphabet_size * sizeof(int));
	// In fact we can reuse bucket_p, but how to compute the
	// l_bucket_end in O(1) space? by +- sign of bucket_p?

	j = new_length;
	for (i = alphabet_size - 1; i > 0; --i) {
		l_bucket_end[i] = j;
		j -= bucket_p[i] - bucket_end[i - 1];
	}
	l_bucket_end[0] = j;

	int *last_deleg =
		(int *)_suffix_array_malloc(alphabet_size * sizeof(int));
	memset(last_deleg, -1, alphabet_size * sizeof(int));

	// init bucket head for filling type S positions
	memcpy(bucket_p, bucket_end, alphabet_size * sizeof(int));

	assert(length > 1);
	// pattern[length] is $, which is the smallest char. The
	// position length - 1 is a type L position.
	c0 = 0;
	b = sorted_suffix + bucket_p[c0];
	int l_c0 = 0;
	int *l_b = sorted_substrings + l_bucket_end[l_c0];
	int s;
	int equal_deleg = 0;
	int new_alphabet_size = new_length;
	for (s = length - 1; s >= 0; --s) {
		i = sorted_suffix[s];
		if (i < 0) {
			i ^= -1;
			if (i >= MASK_N)
				i ^= MASK_N;
			else
				equal_deleg = s;
			--i;
			// update backward
			if ((c1 = pattern[i]) != c0) {
				bucket_p[c0] = b - sorted_suffix;
				b = sorted_suffix + bucket_p[c0 = c1];
			}
			// if i - 1 is type S pos, mark i;
			if (i > 0 && pattern[i - 1] <= c0)
				i ^= -1;

			--b;
			if (equal_deleg == last_deleg[c0])
				*b = i ^ MASK_N;
			else {
				*b = i;
				last_deleg[c0] = equal_deleg;
			}
		}
		else {
			if (i >= MASK_N)
				i ^= MASK_N;
			else
				equal_deleg = s;
			if (i == 0)
				continue;
			// update forward to determine a substring's rank
			--i;
			if ((c1 = pattern[i]) != l_c0) {
				l_bucket_end[l_c0] = l_b - sorted_substrings;
				l_b = sorted_substrings
					+ l_bucket_end[l_c0 = c1];
			}

			// if i - 1 is type S pos, mark i;
			if (i > 0 && pattern[i - 1] < l_c0)
				i ^= -1;

			--l_b;
			if (equal_deleg == last_deleg[l_c0]) {
				*l_b = i ^ MASK_N;
				--new_alphabet_size;
			}
			else {
				*l_b = i;
				last_deleg[l_c0] = equal_deleg;
			}
		}
	}
	i = length - 1;
	if ((c1 = pattern[i]) != l_c0) {
		l_bucket_end[l_c0] = l_b - sorted_substrings;
		l_b = sorted_substrings + l_bucket_end[l_c0 = c1];
	}
	--l_b;
	// if i - 1 is type S pos, mark i;
	if (pattern[i - 1] < l_c0)
		*l_b = ~i;
	else
		*l_b = i;
	free(l_bucket_end);
	free(last_deleg);

	/* compute the new pattern for recursion */
	gettimeofday(&tv, NULL);
	fprintf(stderr, " %d prepare rec at %d.%06d\n",
		length, tv.tv_sec, tv.tv_usec);

	if (new_alphabet_size > 1 && new_alphabet_size < new_length) {
		// save the rank of substrings of each pos in
		// sorted_suffix[pos]. If a position has no tie, its
		// rank ^ -1 will be saved.
		int k = sorted_substrings[new_length - 1];
		if (k < 0)
			k ^= -1;
		//assert(k < MASK_N);
		j = new_alphabet_size - 1;
		for (s = new_length - 2; s >= 0;) {
			i = k;

			do {
				k = sorted_substrings[s];
				if (k < 0)
					k ^= -1;
				if (k >= MASK_N)
					sorted_suffix[k ^ MASK_N] = j;
				else
					break;
				--s;
			} while (s >= 0);
			if (j - 1 == s)
				sorted_suffix[i] = ~j;
			else
				sorted_suffix[i] = j;

			j = s--;
		}
		if (j == 0)
			sorted_suffix[k] = ~j;

		// reduce the new pattern
		int new_r_length = 0;
		for (i = 0; i < new_length; ) {
			s = orig_pos[i];
			if (s < 0)
				k = sorted_suffix[~s];
			else
				k = sorted_suffix[s];
			if (k >= 0) {
				if (s < 0)
					sorted_suffix[~s] = new_r_length;
				else
					sorted_suffix[s] = new_r_length;
				orig_pos[new_r_length++] = s;
				// any char following a char appeared
				// many times in the new pattern
				// should be kept.
				for (++i; (k >= 0) && i < new_length; ++i) {
					s = orig_pos[i];
					if (s < 0) {
						k = sorted_suffix[~s];
						sorted_suffix[~s] = new_r_length;
						// save the position
						// in the new pattern
					}
					else {
						k = sorted_suffix[s];
						sorted_suffix[s] = new_r_length;
					}
					orig_pos[new_r_length++] = s;
				}
			}
			else {
				// mark sorted_suffix[s] to notice that s is removed
				if (s < 0)
					sorted_suffix[~s] = -1;
				else
					sorted_suffix[s] = -1;
				++i;
			}

		}

		// generate the new alphabet and the new pattern.
		int *new_pattern =
			(int *)_suffix_array_malloc(
				new_length * sizeof(int));
			// sorted_suffix; // we can reuse it as new_pattern
			// FIXME. remove test code;
		j = new_alphabet_size + new_r_length - new_length;
		// new_r_alphabet_size
		int *new_bucket_end =
			/*
			(int *)_suffix_array_malloc(
				new_alphabet_size * sizeof(int));
			*/
			sorted_substrings + new_length
			- j; // we can reuse it as new_bucket_end;

		int r_i = new_length;
		// find the first pos which is in the new pattern. it
		// exists.
		for (s = new_length - 1; s >= 0; --s) {
			k = sorted_substrings[s];
			if (k < 0) {
				//assert(k > ~MASK_N);
				i = sorted_suffix[~k];
			}
			else {
				//assert(k < MASK_N);
				i = sorted_suffix[k];
			}
			if (i >= 0)
				break;
			else {
				// this pos is not in the new pattern
				--r_i;
				new_pattern[r_i] = j;
				orig_pos[r_i] = k;
			}
		}
		// assign new_bucket_end[new_r_alphabet_size - 1]
		new_bucket_end[--j] = new_r_length;
		// this pos is in the new pattern
		new_pattern[i] = j;
		--new_r_length;

		for (--s; s >= 0; --s) {
			// the first pos after a in-the-new-pattern pos
			k = sorted_substrings[s];
			if (k < 0)
				if (k <= ~MASK_N) {
					i = sorted_suffix[~k ^ MASK_N];
					//assert(i >= 0);
					goto end;
				}
				else
					i = sorted_suffix[~k];
			else
				if (k >= MASK_N) {
					i = sorted_suffix[k ^ MASK_N];
					//assert(i >= 0);
					goto end;
				}
				else
					i = sorted_suffix[k];
			if (i < 0) {
				--r_i;
				new_pattern[r_i] = j;
				orig_pos[r_i] = k;
				for (--s; s >= 0; --s) {
					k = sorted_substrings[s];
					if (k < 0) {
						//assert(k > ~MASK_N);
						i = sorted_suffix[~k];
					}
					else {
						//assert(k < MASK_N);
						i = sorted_suffix[k];
					}
					if (i >= 0)
						break;
					else {
						// this pos is not in
						// the new pattern
						--r_i;
						new_pattern[r_i] = j;
						orig_pos[r_i] = k;
					}
				}
				if (s < 0)
					break;
			}
			new_bucket_end[--j] = new_r_length;
		end:
			new_pattern[i] = j;
			--new_r_length;
		}
		assert(j == 0);

		*p_new_r_alphabet_size =
			new_alphabet_size + r_i - new_length;
		*p_new_r_length = r_i;
		*p_new_pattern = new_pattern;
		*p_new_bucket_end = new_bucket_end;
		return new_alphabet_size;
	}
	else
		return 0;
}

template <typename T>
void _suffix_array_post_recursion_l(T *pattern, int length,
				    int alphabet_size,
				    int *sorted_suffix,
				    int *bucket_end, int *bucket_p) {
	// init bucket head for filling type S positions
	memcpy(bucket_p, bucket_end, alphabet_size * sizeof(int));

	int i;
	int c0 = 0, c1;
	int *b = sorted_suffix + bucket_p[c0];
	int s;
	for (s = length - 1; s >= 0; --s) {
		i = sorted_suffix[s];
		if (i < 0) {
			i ^= -1;
			sorted_suffix[s] = i;
			--i;
			if ((c1 = pattern[i]) != c0) {
				bucket_p[c0] = b - sorted_suffix;
				b = sorted_suffix + bucket_p[c0 = c1];
			}
			--b;
			// if i - 1 is a type S position, mark i
			if (i > 0 && pattern[i - 1] <= c0)
				*b = ~i;
			else
				*b = i;
		}
	}
}

template <typename T>
void _suffix_array_main_rec(T *pattern, int length,
			    int alphabet_size,
			    int *sorted_suffix,
			    int *bucket_end, int *bucket_p) {
	gettimeofday(&tv, NULL);
	fprintf(stderr, "%d at %d.%06d\n", length, tv.tv_sec, tv.tv_usec);

	// in O2 optimization mode, to decide before saving to
	// orig_pos and bucket is FASTER than to save type S and L
	// positions into orig_pos then decide which one to use and
	// fill in the bucket.

	/* decide whether recur on type L or type S */
	int s_count = 0;

	int same = 0;
	T *c = pattern;
	int c1 = c[0], this_chr;
	int i;
	for (i = 1; i < length; ++i) {
		this_chr = *(++c);
		if (c1 < this_chr) {
			s_count += i - same;
			same = i;
		}
		else if (c1 > this_chr)
			same = i;
		c1 = this_chr;
	}

	int *orig_pos;
	int *sorted_substrings;

	int *new_pattern;
	int new_length;
	int new_alphabet_size;
	int *new_bucket_end;
	int *new_bucket_p;
	int *new_sorted_suffix;

	int s;
	int *b;

	int new_r_length;
	int new_r_alphabet_size = 0;
	if (s_count < ((length + 1) >> 1)) {
		fprintf(stderr, " %d S is smaller\n", length);
		new_length = s_count;
		orig_pos =
			(int *)_suffix_array_malloc(new_length * sizeof(int));
		sorted_substrings =
			(int *)_suffix_array_malloc(new_length * sizeof(int));
		new_alphabet_size =
			_suffix_array_sort_substring_s<T>(
				pattern, length, new_length,
				alphabet_size,
				orig_pos, sorted_suffix,
				bucket_end, bucket_p,
				sorted_substrings,
				&new_pattern, &new_bucket_end,
				&new_r_alphabet_size, &new_r_length);
		// init bucket for filling the type S positions.
		memcpy(bucket_p, bucket_end, alphabet_size * sizeof(int));
		c1 = alphabet_size - 1;
		b = sorted_suffix + bucket_p[c1];
		// check if recursion and run the recursion
		if (new_alphabet_size > 0) {
			new_sorted_suffix = sorted_suffix;
			new_bucket_p =
				(int *)_suffix_array_malloc(
					new_r_alphabet_size * sizeof(int));
			_suffix_array_main_rec<int>(
				new_pattern, new_r_length,
				new_r_alphabet_size,
				new_sorted_suffix,
				new_bucket_end, new_bucket_p);
			free(new_bucket_p);

			gettimeofday(&tv, NULL);
			fprintf(stderr, " %d post rec start at %d.%06d\n",
				length, tv.tv_sec, tv.tv_usec);

			// no need to free new_pattern and new_bucket_end

			// WARNING: we can only reuse the
			// sorted_suffix from the beginning as
			// new_sorted_suffix

			// merge reduced pos and the original
			// positions of the chars in the new_pattern.
			int t = new_length - 1;
			s = new_r_length - 1;
			int tmp;

			int j;
			while (t >= new_r_length) {
				for (j = new_pattern[t]; // must be valid
				     s >= 0 && new_pattern[
					     tmp = new_sorted_suffix[s]] >= j;
				     --s) {
					i = orig_pos[tmp];
					if (i < 0)
						this_chr = pattern[~i];
					else
						this_chr = pattern[i];
					if (this_chr != c1) {
						bucket_p[c1] = b - sorted_suffix;
						b = sorted_suffix
							+ bucket_p[c1 = this_chr];
					}
					*(--b) = i;
				}
				if (s < 0)
					break;
				for (j = new_pattern[tmp];
				     t >= new_r_length && new_pattern[t] > j;
				     --t) {
					i = orig_pos[t];
					if (i < 0)
						this_chr = pattern[~i];
					else
						this_chr = pattern[i];
					if (this_chr != c1) {
						bucket_p[c1] = b - sorted_suffix;
						b = sorted_suffix
							+ bucket_p[c1 = this_chr];
					}
					*(--b) = i;
				}
			}
			for (; t >= new_r_length; --t) {
				i = orig_pos[t];
				if (i < 0)
					this_chr = pattern[~i];
				else
					this_chr = pattern[i];
				if (this_chr != c1) {
					bucket_p[c1] = b - sorted_suffix;
					b = sorted_suffix
						+ bucket_p[c1 = this_chr];
				}
				*(--b) = i;
			}
			for (; s >= 0; --s) {
				i = orig_pos[new_sorted_suffix[s]];
				if (i < 0)
					this_chr = pattern[~i];
				else
					this_chr = pattern[i];
				if (this_chr != c1) {
					bucket_p[c1] = b - sorted_suffix;
					b = sorted_suffix
						+ bucket_p[c1 = this_chr];
				}
				*(--b) = i;
			}

			free(new_pattern);
			// FIXME. remove test code
		}
		else {
			gettimeofday(&tv, NULL);
			fprintf(stderr, " %d post rec start at %d.%06d\n",
				length, tv.tv_sec, tv.tv_usec);

			for (s = new_length - 1; s >= 0; --s) {
				i = sorted_substrings[s];
				if (i < 0) {
					if (i <= ~MASK_N)
						i ^= MASK_N;
					this_chr = pattern[~i];
				}
				else {
					if (i >= MASK_N)
						i ^= MASK_N;
					this_chr = pattern[i];
				}
				if (this_chr != c1) {
					bucket_p[c1] = b - sorted_suffix;
					b = sorted_suffix
						+ bucket_p[c1 = this_chr];
				}
				*(--b) = i;
			}
		}
		free(sorted_substrings);
		free(orig_pos);

		gettimeofday(&tv, NULL);
		fprintf(stderr, " %d post rec finish at %d.%06d\n",
			length, tv.tv_sec, tv.tv_usec);

		_suffix_array_post_recursion_s<T>(pattern, length,
						  alphabet_size,
						  sorted_suffix,
						  bucket_end, bucket_p);
	}
	else {
		fprintf(stderr, " %d L is smaller\n", length);
		new_length = length - s_count;
		orig_pos =
			(int *)_suffix_array_malloc(new_length * sizeof(int));
		sorted_substrings =
			(int *)_suffix_array_malloc(new_length * sizeof(int));
		new_alphabet_size =
			_suffix_array_sort_substring_l<T>(
				pattern, length, new_length,
				alphabet_size,
				orig_pos, sorted_suffix,
				bucket_end, bucket_p,
				sorted_substrings,
				&new_pattern, &new_bucket_end,
				&new_r_alphabet_size, &new_r_length);
		// init bucket for filling the type L positions.
		bucket_p[0] = 0;
		memcpy(bucket_p + 1, bucket_end,
		       (alphabet_size - 1) * sizeof(int));
		c1 = 0;
		b = sorted_suffix + bucket_p[c1];
		// check if recursion and run the recursion
		if (new_alphabet_size > 0) {
			new_sorted_suffix = sorted_suffix
				+ length - new_r_length;
			new_bucket_p =
				(int *)_suffix_array_malloc(
					new_r_alphabet_size * sizeof(int));
			_suffix_array_main_rec<int>(
				new_pattern, new_r_length,
				new_r_alphabet_size,
				new_sorted_suffix,
				new_bucket_end, new_bucket_p);
			free(new_bucket_p);

			gettimeofday(&tv, NULL);
			fprintf(stderr, " %d post rec start at %d.%06d\n",
				length, tv.tv_sec, tv.tv_usec);

			// no need to free new_pattern and new_bucket_end

			// WARNING: we can only reuse the end side of
			// sorted_suffix as new_sorted_suffix

			// merge reduced pos and the original
			// positions of the chars in the new_pattern.
			int t = new_r_length;
			s = 0;
			int tmp;

			int j;
			while (t < new_length) {
				for (j = new_pattern[t];
				     s < new_r_length && new_pattern[
					     tmp = new_sorted_suffix[s]] < j;
				     ++s) {
					i = orig_pos[tmp];
					if (i < 0)
						this_chr = pattern[~i];
					else
						this_chr = pattern[i];
					if (this_chr != c1) {
						bucket_p[c1] = b - sorted_suffix;
						b = sorted_suffix
							+ bucket_p[c1 = this_chr];
					}
					*(b++) = i;
				}
				if (s == new_r_length)
					break;
				for (j = new_pattern[tmp];
				     t < new_length && new_pattern[t] <= j;
				     ++t) {
					i = orig_pos[t];
					if (i < 0)
						this_chr = pattern[~i];
					else
						this_chr = pattern[i];
					if (this_chr != c1) {
						bucket_p[c1] = b - sorted_suffix;
						b = sorted_suffix
							+ bucket_p[c1 = this_chr];
					}
					*(b++) = i;
				}
			}
			for (; t < new_length; ++t) {
				i = orig_pos[t];
				if (i < 0)
					this_chr = pattern[~i];
				else
					this_chr = pattern[i];
				if (this_chr != c1) {
					bucket_p[c1] = b - sorted_suffix;
					b = sorted_suffix
						+ bucket_p[c1 = this_chr];
				}
				*(b++) = i;
			}
			for (; s < new_r_length; ++s) {
				i = orig_pos[new_sorted_suffix[s]];
				if (i < 0)
					this_chr = pattern[~i];
				else
					this_chr = pattern[i];
				if (this_chr != c1) {
					bucket_p[c1] = b - sorted_suffix;
					b = sorted_suffix
						+ bucket_p[c1 = this_chr];
				}
				*(b++) = i;
			}

			free(new_pattern);
			// FIXME. remove test code
		}
		else {
			gettimeofday(&tv, NULL);
			fprintf(stderr, " %d post rec start at %d.%06d\n",
				length, tv.tv_sec, tv.tv_usec);

			for (s = 0; s < new_length; ++s) {
				i = sorted_substrings[s];
				if (i < 0) {
					if (i <= ~MASK_N)
						i ^= MASK_N;
					this_chr = pattern[~i];
				}
				else {
					if (i >= MASK_N)
						i ^= MASK_N;
					this_chr = pattern[i];
				}
				if (this_chr != c1) {
					bucket_p[c1] = b - sorted_suffix;
					b = sorted_suffix
						+ bucket_p[c1 = this_chr];
				}
				*(b++) = i;
			}
		}
		free(sorted_substrings);
		free(orig_pos);

		gettimeofday(&tv, NULL);
		fprintf(stderr, " %d post rec finish at %d.%06d\n",
			length, tv.tv_sec, tv.tv_usec);

		_suffix_array_post_recursion_l<T>(pattern, length,
						  alphabet_size,
						  sorted_suffix,
						  bucket_end, bucket_p);
	}

	gettimeofday(&tv, NULL);
	fprintf(stderr, "end of %d at %d.%06d\n", length, tv.tv_sec, tv.tv_usec);
}

// this is the no reduction version. Reduce the new_pattern is useless
// in the first level.
template <typename T>
int _suffix_array_sort_substring_s_no_reduction(
	T *pattern, int length, int new_length,	int alphabet_size,
	int *orig_pos, int *sorted_suffix, int *bucket_end, int *bucket_p,
	int *sorted_substrings, int **_new_pattern, int **_new_bucket_end) {
	/* fill the type S positions into the bucket */
	memcpy(bucket_p, bucket_end, alphabet_size * sizeof(int));
	//memset(sorted_suffix, 0, length * sizeof(int));
	int i;
	int j = 0;
	int same = 0; // same means the same char
	T *c = pattern;
	int c1 = c[0], c0 = c1, this_chr; // we should avoid overflow of signed char
	int *b = sorted_suffix + bucket_p[c0];
	int r = 0;
	// assume -1 is a type S position. That means, if 0 is a type
	// S position, it won't enter the forward update branch.
	for (i = 1; i < length; ++i) {
		this_chr = *(++c);
		if (c1 < this_chr) {
			if (c1 != c0) {
				bucket_p[c0] = b - sorted_suffix;
				b = sorted_suffix + bucket_p[c0 = c1];
			}

			--b;
			if (r == 1) {
				// same - 1 is a type L position, mark same.
				orig_pos[j] = ~same;
				*b = ~same ^ MASK_N;
				r = 0;
			}
			else {
				orig_pos[j] = same;
				*b = same ^ MASK_N;
			}
			++j;

			for (++same; same < i; ++same) {
				orig_pos[j++] = same;
				*(--b) = same ^ MASK_N;
			}
		}
		else if (c1 > this_chr) {
			same = i;
			r = 1;
		}
		c1 = this_chr;
	}
	// if same != length - 1, then the last same block is consist
	// of all type L positions.
	bucket_p[c0] = b - sorted_suffix;
	// alphabet_size is at least 2
	for (i = alphabet_size - 1; i >= 0; --i)
		if ((j = bucket_p[i]) < bucket_end[i])
			sorted_suffix[j] ^= MASK_N;
	/* no need to assign equal value for type L suffixes and no
	 * need to set the order for type L suffixes. */

	/* sort the substrings */
	int *s_bucket_begin =
		(int *)_suffix_array_malloc(alphabet_size * sizeof(int));
	// In fact we can reuse bucket_p, but how to compute the
	// s_bucket_begin in O(1) space? by +- sign of bucket_p?

	j = 0;
	for (i = 0; i < alphabet_size; ++i) {
		s_bucket_begin[i] = j;
		j += bucket_end[i] - bucket_p[i];
	}

	int *last_deleg =
		(int *)_suffix_array_malloc(alphabet_size * sizeof(int));
	memset(last_deleg, -1, alphabet_size * sizeof(int));

	// init bucket head for filling type L positions
	bucket_p[0] = 0;
	memcpy(bucket_p + 1, bucket_end, (alphabet_size - 1) * sizeof(int));

	assert(length > 1);
	// pattern[length] is $, which is the smallest char. The
	// position length - 1 is a type L position.
	i = length - 1;
	c0 = pattern[i];
	b = sorted_suffix + bucket_p[c0];
	// if i - 1 is type L pos, mark i;
	if (pattern[i - 1] >= c0)
		*b = ~i;
	else
		*b = i;
	++b;
	int s_c0 = 0;
	int *s_b = sorted_substrings + s_bucket_begin[s_c0];
	int s;
	int equal_deleg = 0;
	int new_alphabet_size = new_length;
	for (s = 0; s < length; ++s) {
		i = sorted_suffix[s];
		if (i < 0) {
			i ^= -1;
			if (i >= MASK_N)
				i ^= MASK_N;
			else
				equal_deleg = s;
			--i;
			// update forward
			if ((c1 = pattern[i]) != c0) {
				bucket_p[c0] = b - sorted_suffix;
				b = sorted_suffix + bucket_p[c0 = c1];
			}
			// if i - 1 is type L pos, mark i;
			if (i > 0 && pattern[i - 1] >= c0)
				i ^= -1;

			if (equal_deleg == last_deleg[c0])
				*b = i ^ MASK_N;
			else {
				*b = i;
				last_deleg[c0] = equal_deleg;
			}
			++b;
		}
		else {
			if (i >= MASK_N)
				i ^= MASK_N;
			else
				equal_deleg = s;
			if (i == 0)
				continue;
			// update backward to determine a substring's rank
			--i;
			if ((c1 = pattern[i]) != s_c0) {
				s_bucket_begin[s_c0] = s_b - sorted_substrings;
				s_b = sorted_substrings
					+ s_bucket_begin[s_c0 = c1];
			}

			// if i - 1 is type L pos, mark i;
			if (i > 0 && pattern[i - 1] > s_c0)
				i ^= -1;

			if (equal_deleg == last_deleg[s_c0]) {
				*s_b = i ^ MASK_N;
				--new_alphabet_size;
			}
			else {
				*s_b = i;
				last_deleg[s_c0] = equal_deleg;
			}
			++s_b;
		}
	}
	free(s_bucket_begin);
	free(last_deleg);

	/* compute the new pattern for recursion */
	gettimeofday(&tv, NULL);
	fprintf(stderr, " %d prepare rec at %d.%06d\n",
		length, tv.tv_sec, tv.tv_usec);

	if (new_alphabet_size > 1 && new_alphabet_size < new_length) {
		// no reduction
		int *new_bucket_end =
			/*
			(int *)_suffix_array_malloc(
				new_alphabet_size * sizeof(int));
			*/
			sorted_substrings; // we can reuse it as new_bucket_end;

		i = sorted_substrings[0];
		if (i < 0)
			i ^= -1;
		assert(i < MASK_N);
		j = 0;
		sorted_suffix[i] = j;
		for (s = 1; s < new_length; ++s) {
			i = sorted_substrings[s];
			if (i < 0)
				i ^= -1;
			if (i < MASK_N) {
				new_bucket_end[j] = s;
				++j;
			}
			else
				i ^= MASK_N;
			sorted_suffix[i] = j;
		}
		new_bucket_end[j] = s;
		int *new_pattern =
			(int *)_suffix_array_malloc(
				new_length * sizeof(int));
			// sorted_suffix; // we can reuse it as new_pattern
			// FIXME. remove test code;
		for (i = 0; i < new_length; ++i) {
			s = orig_pos[i];
			if (s < 0)
				new_pattern[i] = sorted_suffix[~s];
			else
				new_pattern[i] = sorted_suffix[s];
		}

		*_new_pattern = new_pattern;
		*_new_bucket_end = new_bucket_end;
		return new_alphabet_size;
	}
	else
		return 0;
}

template <typename T>
int _suffix_array_sort_substring_l_no_reduction(
	T *pattern, int length, int new_length, int alphabet_size,
	int *orig_pos, int *sorted_suffix, int *bucket_end, int *bucket_p,
	int *sorted_substrings, int **_new_pattern, int **_new_bucket_end) {
	/* fill the type L positions into the bucket */
	bucket_p[0] = 0;
	memcpy(bucket_p + 1, bucket_end, (alphabet_size - 1) * sizeof(int));
	//memset(sorted_suffix, 0, length * sizeof(int));
	int i;
	int j = new_length;
	T *c = pattern + length - 1;
	int c1 = c[0], c0 = c1, this_chr;
	int *b = sorted_suffix + bucket_p[c0];
	int r = 1;
	// assume -1 is a type L position. That means, if 0 is a type
	// L position, it won't enter the backward update branch.
	for (i = length - 1; i > 0; --i) {
		this_chr = *(--c);
		if (this_chr > c1) {
			// this_chr is L
			if (r == 1) {
				// the next chr is L, save it
				orig_pos[--j] = i;
				*(b++) = i ^ MASK_N;
				bucket_p[c0] = b - sorted_suffix;
				b = sorted_suffix + bucket_p[c0 = this_chr];
			}
			else {
				if (c1 != c0) {
					bucket_p[c0] = b - sorted_suffix;
					b = sorted_suffix
						+ bucket_p[c0 = this_chr];
				}
				r = 1;
			}
		}
		else if (this_chr < c1) {
			// this_chr is S
			if (r == 1) {
				// the next chr is L, save it, with mark
				orig_pos[--j] = ~i;
				*(b++) = ~i ^ MASK_N;
				r = 0;
			}
		}
		else {
			if (r == 1) {
				// the next chr is L, this chr is also
				// L, save the next chr
				orig_pos[--j] = i;
				*(b++) = i ^ MASK_N;
			}
		}
		c1 = this_chr;
	}
	// if the first chr is L, save it.
	if (r == 1) {
		orig_pos[--j] = i;
		*(b++) = i ^ MASK_N;
	}
	bucket_p[c0] = b - sorted_suffix;
	// alphabet_size is at least 1
	if ((j = bucket_p[0] - 1) >= 0)
		sorted_suffix[j] ^= MASK_N;
	for (i = 1; i < alphabet_size; ++i)
		if ((j = bucket_p[i] - 1) >= bucket_end[i - 1])
			sorted_suffix[j] ^= MASK_N;
	/* no need to assign equal value for type L suffixes and no
	 * need to set the order for type L suffixes. */

	/* sort the substrings */
	int *l_bucket_end =
		(int *)_suffix_array_malloc(alphabet_size * sizeof(int));
	// In fact we can reuse bucket_p, but how to compute the
	// l_bucket_end in O(1) space? by +- sign of bucket_p?

	j = new_length;
	for (i = alphabet_size - 1; i > 0; --i) {
		l_bucket_end[i] = j;
		j -= bucket_p[i] - bucket_end[i - 1];
	}
	l_bucket_end[0] = j;

	int *last_deleg =
		(int *)_suffix_array_malloc(alphabet_size * sizeof(int));
	memset(last_deleg, -1, alphabet_size * sizeof(int));

	// init bucket head for filling type S positions
	memcpy(bucket_p, bucket_end, alphabet_size * sizeof(int));

	assert(length > 1);
	// pattern[length] is $, which is the smallest char. The
	// position length - 1 is a type L position.
	c0 = 0;
	b = sorted_suffix + bucket_p[c0];
	int l_c0 = 0;
	int *l_b = sorted_substrings + l_bucket_end[l_c0];
	int s;
	int equal_deleg = 0;
	int new_alphabet_size = new_length;
	for (s = length - 1; s >= 0; --s) {
		i = sorted_suffix[s];
		if (i < 0) {
			i ^= -1;
			if (i >= MASK_N)
				i ^= MASK_N;
			else
				equal_deleg = s;
			--i;
			// update backward
			if ((c1 = pattern[i]) != c0) {
				bucket_p[c0] = b - sorted_suffix;
				b = sorted_suffix + bucket_p[c0 = c1];
			}
			// if i - 1 is type S pos, mark i;
			if (i > 0 && pattern[i - 1] <= c0)
				i ^= -1;

			--b;
			if (equal_deleg == last_deleg[c0])
				*b = i ^ MASK_N;
			else {
				*b = i;
				last_deleg[c0] = equal_deleg;
			}
		}
		else {
			if (i >= MASK_N)
				i ^= MASK_N;
			else
				equal_deleg = s;
			if (i == 0)
				continue;
			// update forward to determine a substring's rank
			--i;
			if ((c1 = pattern[i]) != l_c0) {
				l_bucket_end[l_c0] = l_b - sorted_substrings;
				l_b = sorted_substrings
					+ l_bucket_end[l_c0 = c1];
			}

			// if i - 1 is type S pos, mark i;
			if (i > 0 && pattern[i - 1] < l_c0)
				i ^= -1;

			--l_b;
			if (equal_deleg == last_deleg[l_c0]) {
				*l_b = i ^ MASK_N;
				--new_alphabet_size;
			}
			else {
				*l_b = i;
				last_deleg[l_c0] = equal_deleg;
			}
		}
	}
	i = length - 1;
	if ((c1 = pattern[i]) != l_c0) {
		l_bucket_end[l_c0] = l_b - sorted_substrings;
		l_b = sorted_substrings + l_bucket_end[l_c0 = c1];
	}
	--l_b;
	// if i - 1 is type S pos, mark i;
	if (pattern[i - 1] < l_c0)
		*l_b = ~i;
	else
		*l_b = i;
	free(l_bucket_end);
	free(last_deleg);

	/* compute the new pattern for recursion */
	gettimeofday(&tv, NULL);
	fprintf(stderr, " %d prepare rec at %d.%06d\n",
		length, tv.tv_sec, tv.tv_usec);

	if (new_alphabet_size > 1 && new_alphabet_size < new_length) {
		// no reduction
		int *new_bucket_end =
			/*
			(int *)_suffix_array_malloc(
				new_alphabet_size * sizeof(int));
			*/
			sorted_substrings + new_length
			- new_alphabet_size; // we can reuse it as
					     // new_bucket_end;

		i = sorted_substrings[new_length - 1];
		if (i < 0)
			i ^= -1;
		assert(i < MASK_N);
		j = new_alphabet_size - 1;
		new_bucket_end[j] = new_length;
		sorted_suffix[i] = j;
		for (s = new_length - 2; s >= 0; --s) {
			i = sorted_substrings[s];
			if (i < 0)
				i ^= -1;
			if (i < MASK_N)
				new_bucket_end[--j] = s + 1;
			else
				i ^= MASK_N;
			sorted_suffix[i] = j;
		}
		int *new_pattern =
			(int *)_suffix_array_malloc(
				new_length * sizeof(int));
			// sorted_suffix; // we can reuse it as new_pattern
			// FIXME. remove test code;
		for (i = 0; i < new_length; ++i) {
			s = orig_pos[i];
			if (s < 0)
				new_pattern[i] = sorted_suffix[~s];
			else
				new_pattern[i] = sorted_suffix[s];
		}

		*_new_pattern = new_pattern;
		*_new_bucket_end = new_bucket_end;
		return new_alphabet_size;
	}
	else
		return 0;
}

template <typename T>
void _suffix_array_main(T *pattern, int length,
			int alphabet_size,
			int *sorted_suffix,
			int *bucket_end, int *bucket_p) {
	gettimeofday(&tv, NULL);
	fprintf(stderr, "%d at %d.%06d\n", length, tv.tv_sec, tv.tv_usec);

	// in O2 optimization mode, to decide before saving to
	// orig_pos and bucket is FASTER than to save type S and L
	// positions into orig_pos then decide which one to use and
	// fill in the bucket.

	/* decide whether recur on type L or type S */
	int s_count = 0;

	int same = 0;
	T *c = pattern;
	int c1 = c[0], this_chr;
	int i;
	for (i = 1; i < length; ++i) {
		this_chr = *(++c);
		if (c1 < this_chr) {
			s_count += i - same;
			same = i;
		}
		else if (c1 > this_chr)
			same = i;
		c1 = this_chr;
	}

	int *orig_pos;
	int *sorted_substrings;

	int *new_pattern;
	int new_length;
	int *new_bucket_end;
	int *new_bucket_p;
	int *new_sorted_suffix;

	int s;
	int *b;
	if (s_count < ((length + 1) >> 1)) {
		fprintf(stderr, " %d S is smaller\n", length);
		new_length = s_count;
		orig_pos =
			(int *)_suffix_array_malloc(new_length * sizeof(int));
		sorted_substrings =
			(int *)_suffix_array_malloc(new_length * sizeof(int));
		int new_alphabet_size =
			_suffix_array_sort_substring_s_no_reduction<T>(
				pattern, length, new_length,
				alphabet_size,
				orig_pos, sorted_suffix,
				bucket_end, bucket_p,
				sorted_substrings,
				&new_pattern, &new_bucket_end);
		// init bucket for filling the type S positions.
		memcpy(bucket_p, bucket_end, alphabet_size * sizeof(int));
		c1 = alphabet_size - 1;
		b = sorted_suffix + bucket_p[c1];
		// check if recursion and run the recursion
		if (new_alphabet_size > 0) {
			new_sorted_suffix = sorted_suffix;
			new_bucket_p =
				(int *)_suffix_array_malloc(
					new_alphabet_size * sizeof(int));
			_suffix_array_main_rec<int>(
				new_pattern, new_length,
				new_alphabet_size,
				new_sorted_suffix,
				new_bucket_end, new_bucket_p);
			free(new_bucket_p);

			free(new_pattern);
			// FIXME. remove test code

			gettimeofday(&tv, NULL);
			fprintf(stderr, " %d post rec start at %d.%06d\n",
				length, tv.tv_sec, tv.tv_usec);

			// no need to free new_pattern and new_bucket_end

			// WARNING: we can only reuse the
			// sorted_suffix from the beginning as
			// new_sorted_suffix
			for (s = new_length - 1; s >= 0; --s) {
				i = orig_pos[new_sorted_suffix[s]];
				if (i < 0)
					this_chr = pattern[~i];
				else
					this_chr = pattern[i];
				if (this_chr != c1) {
					bucket_p[c1] = b - sorted_suffix;
					b = sorted_suffix
						+ bucket_p[c1 = this_chr];
				}
				*(--b) = i;
			}
		}
		else {
			gettimeofday(&tv, NULL);
			fprintf(stderr, " %d post rec start at %d.%06d\n",
				length, tv.tv_sec, tv.tv_usec);

			for (s = new_length - 1; s >= 0; --s) {
				i = sorted_substrings[s];
				if (i < 0) {
					if (i & MASK_N == 0)
						i ^= MASK_N;
					this_chr = pattern[~i];
				}
				else {
					if (i >= MASK_N)
						i ^= MASK_N;
					this_chr = pattern[i];
				}
				if (this_chr != c1) {
					bucket_p[c1] = b - sorted_suffix;
					b = sorted_suffix
						+ bucket_p[c1 = this_chr];
				}
				*(--b) = i;
			}
		}
		free(sorted_substrings);
		free(orig_pos);

		gettimeofday(&tv, NULL);
		fprintf(stderr, " %d post rec finish at %d.%06d\n",
			length, tv.tv_sec, tv.tv_usec);

		_suffix_array_post_recursion_s<T>(pattern, length,
						  alphabet_size,
						  sorted_suffix,
						  bucket_end, bucket_p);
	}
	else {
		fprintf(stderr, " %d L is smaller\n", length);
		new_length = length - s_count;
		orig_pos =
			(int *)_suffix_array_malloc(new_length * sizeof(int));
		sorted_substrings =
			(int *)_suffix_array_malloc(new_length * sizeof(int));
		int new_alphabet_size =
			_suffix_array_sort_substring_l_no_reduction<T>(
				pattern, length, new_length,
				alphabet_size,
				orig_pos, sorted_suffix,
				bucket_end, bucket_p,
				sorted_substrings,
				&new_pattern, &new_bucket_end);
		// init bucket for filling the type L positions.
		bucket_p[0] = 0;
		memcpy(bucket_p + 1, bucket_end,
		       (alphabet_size - 1) * sizeof(int));
		c1 = 0;
		b = sorted_suffix + bucket_p[c1];
		// check if recursion and run the recursion
		if (new_alphabet_size > 0) {
			new_sorted_suffix = sorted_suffix
				+ length - new_length;
			new_bucket_p =
				(int *)_suffix_array_malloc(
					new_alphabet_size * sizeof(int));
			_suffix_array_main_rec<int>(
				new_pattern, new_length,
				new_alphabet_size,
				new_sorted_suffix,
				new_bucket_end, new_bucket_p);
			free(new_bucket_p);

			free(new_pattern);
			// FIXME. remove test code

			gettimeofday(&tv, NULL);
			fprintf(stderr, " %d post rec start at %d.%06d\n",
				length, tv.tv_sec, tv.tv_usec);

			// no need to free new_pattern and new_bucket_end

			// WARNING: we can only reuse the end side of
			// sorted_suffix as new_sorted_suffix
			for (s = 0; s < new_length; ++s) {
				i = orig_pos[new_sorted_suffix[s]];
				if (i < 0)
					this_chr = pattern[~i];
				else
					this_chr = pattern[i];
				if (this_chr != c1) {
					bucket_p[c1] = b - sorted_suffix;
					b = sorted_suffix
						+ bucket_p[c1 = this_chr];
				}
				*(b++) = i;
			}
		}
		else {
			gettimeofday(&tv, NULL);
			fprintf(stderr, " %d post rec start at %d.%06d\n",
				length, tv.tv_sec, tv.tv_usec);

			for (s = 0; s < new_length; ++s) {
				i = sorted_substrings[s];
				if (i < 0) {
					if (i & MASK_N == 0)
						i ^= MASK_N;
					this_chr = pattern[~i];
				}
				else {
					if (i >= MASK_N)
						i ^= MASK_N;
					this_chr = pattern[i];
				}
				if (this_chr != c1) {
					bucket_p[c1] = b - sorted_suffix;
					b = sorted_suffix
						+ bucket_p[c1 = this_chr];
				}
				*(b++) = i;
			}
		}
		free(sorted_substrings);
		free(orig_pos);

		gettimeofday(&tv, NULL);
		fprintf(stderr, " %d post rec finish at %d.%06d\n",
			length, tv.tv_sec, tv.tv_usec);

		_suffix_array_post_recursion_l<T>(pattern, length,
						  alphabet_size,
						  sorted_suffix,
						  bucket_end, bucket_p);
	}

	gettimeofday(&tv, NULL);
	fprintf(stderr, "end of %d at %d.%06d\n", length, tv.tv_sec, tv.tv_usec);
}

void suffix_array(char *pattern, int length) {
	int alphabet_size = 256;
	int *sorted_suffix = (int *)_suffix_array_malloc(length * sizeof(int));
	int *bucket_end = (int *)_suffix_array_malloc(alphabet_size * sizeof(int));
	int *bucket_p = (int *)_suffix_array_malloc(alphabet_size * sizeof(int));

	assert(length > 2);

	_suffix_array_init((unsigned char *)pattern, length, bucket_end);
	_suffix_array_main<unsigned char>((unsigned char *)pattern, length, alphabet_size,
					  sorted_suffix,
					  bucket_end, bucket_p);
	free(bucket_end);
	free(bucket_p);
	int i;
	for (i = 0; i < length; ++i)
		printf("%d\n", sorted_suffix[i]);
	free(sorted_suffix);
}

#define MAXN 30000002

char s[MAXN];

int main () {
	fgets(s, MAXN, stdin);
	int len = strlen(s);
	if (s[len - 2] == '\r')
		--len;
	s[--len] = 0;
	suffix_array(s, len);
	return 0;
}
