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
int _suffix_array_init(char *pattern, int length,
		       int alphabet_size,
		       int *sorted_suffix,
		       int *bucket_end, int *bucket_p,
		       int **new_orig_pos) {
	gettimeofday(&tv, NULL);
	fprintf(stderr, "init %d at %d.%06d\n", length, tv.tv_sec, tv.tv_usec);

	/* construct the bucket */
	// char_count must be inited to 0
	int *char_count = (int *)_suffix_array_malloc(alphabet_size * sizeof(int));
	memset(char_count, 0, alphabet_size * sizeof(int));
	char *p = pattern + length;
	char *c;
	for (c = pattern; c < p; ++c)
		++char_count[*c];
	int tmp_int1 = 0;
	int i;
	for (i = 0; i < alphabet_size; ++i) {
		/* tmp_int1 : bucket_end[i] */
		tmp_int1 += char_count[i];
		bucket_end[i] = tmp_int1;
	}

	memset(sorted_suffix, 0, length * sizeof(int));

	/* find all type S* positions and put them into their
	 * bucket */
	int *_new_orig_pos =
		(int *)_suffix_array_malloc(length * sizeof(int)); // FIXME
	*new_orig_pos = _new_orig_pos;

	int ret = 0;

	memcpy(bucket_p, bucket_end, alphabet_size * sizeof(int));
	int same = 0;
	c = pattern;
	char c1 = c[0], this_chr, c0 = c1;
	int *b = sorted_suffix + bucket_p[c0];
	char r = 0;
	for (i = 1; i < length; ++i) {
		this_chr = *(++c);
		if (c1 < this_chr) {
			// type[same:i - 1] is S
			if (r == 1) {
				// same is a type S* position
				_new_orig_pos[ret++] = same;
				if (c1 != c0) {
					bucket_p[c0] = b - sorted_suffix;
					b = sorted_suffix
						+ bucket_p[c0 = c1];
				}
				*(--b) = ~same ^ MASK_N;
				// if it's equal to its next, xor the MASK_N bit.
				r = 0;
			}
			//same = i; // can be omitted
		}
		else if (c1 > this_chr) {
			// type[same:i - 1] is L
			r = 1;
			same = i;
		}
		c1 = this_chr;
	}
	bucket_p[c0] = b - sorted_suffix;

	// alphabet_size is at least 2
	for (i = alphabet_size - 1; i >= 0; --i)
		if ((tmp_int1 = bucket_p[i]) < bucket_end[i])
			sorted_suffix[tmp_int1] ^= MASK_N;
	/* no need to assign equal value for type L suffixes and no
	 * need to set the order for type L suffixes. */

	free(char_count);
	return ret;
}

template <typename T>
void _suffix_array_worker(T *pattern, int length, int new_length,
			  int *new_orig_pos,
			  int alphabet_size,
			  int *sorted_suffix,
			  int *bucket_end, int *bucket_p) {
	gettimeofday(&tv, NULL);
	fprintf(stderr, "%d at %d.%06d\n", length, tv.tv_sec, tv.tv_usec);

	/* sort substrings */
	int s;
	int tmp_int1;
	int tmp_int2;
	int i;
	int c0, c1;
	int *b;
	// first round
	bucket_p[0] = 0;
	memcpy(bucket_p + 1, bucket_end, (alphabet_size - 1) * sizeof(int));

	i = length - 1;
	c0 = pattern[i];
	b = sorted_suffix + bucket_p[c0];
	if (pattern[i - 1] < c0)
		*b = i;
	else
		*b = ~i;
	++b;

	int *last_deleg =
		(int *)_suffix_array_malloc(alphabet_size * sizeof(int));
	memset(last_deleg, -1, alphabet_size * sizeof(int));
	int equal_deleg = 0;
	for (s = 0; s < length; ++s) {
		i = sorted_suffix[s];
		if (i < 0) {
			i = ~i;
		 	if (i >= MASK_N) {
				i ^= MASK_N;
				sorted_suffix[s] = MASK_N;
			}
			else {
				equal_deleg = s;
				sorted_suffix[s] = 0;
			}
			// The value for type LL/LS has been
			// erased. Then in the second round, we can
			// identify type S* (LS) position by checking
			// if its value is < 0.
			// FIXME: consider defer the above MASK_N?
			--i;
			if ((c1 = pattern[i]) != c0) {
				bucket_p[c0] = b - sorted_suffix;
				b = sorted_suffix + bucket_p[c0 = c1];
			}
			if (i > 0 && pattern[i - 1] >= c0)
				i ^= -1;
			// set 0 (if 0 is a type L pos) as 0.
			if (equal_deleg == last_deleg[c0])
				*b = i ^ MASK_N;
 			else {
				*b = i;
				last_deleg[c0] = equal_deleg;
			}
			++b;
		}
		else if (i < MASK_N)
			equal_deleg = s;
	}

	// second round
	/* init new_ series variables for the recursion */
	int new_alphabet_size = new_length;
	int *new_pattern = (int *)_suffix_array_malloc(new_length * sizeof(int));
	int *new_sorted_suffix =
		(int *)_suffix_array_malloc(new_length * sizeof(int));
	int *new_rank_of = (int *)_suffix_array_malloc(new_length * sizeof(int));
	int tmp_int3;

 	memcpy(bucket_p, bucket_end, alphabet_size * sizeof(int));

	memset(last_deleg, 0, alphabet_size * sizeof(int));
	b = sorted_suffix + bucket_p[c0];
	equal_deleg = length - 1;
	tmp_int2 = 1;
	int j = new_length;
	int last_sstar_deleg = length; // obviously != equal_deleg
	for (s = length - 1; s >= 0; --s) {
		i = sorted_suffix[s];
		// how to recognize type S* pos?
		if (i >= 0) {
			if (i >= MASK_N) {
				i ^= MASK_N;
				tmp_int2 = 0;
			}
			if (i > 0) {
				--i;
				if ((c1 = pattern[i]) != c0) {
					bucket_p[c0] = b - sorted_suffix;
					b = sorted_suffix + bucket_p[c0 = c1];
				}
				if (equal_deleg == last_deleg[c0])
					*b ^= MASK_N;
				else
					last_deleg[c0] = equal_deleg;
				--b;
				if (i == 0 || pattern[i - 1] <= c0)
					*b = i;
				else
					*b = ~i;
			}
		}
		else {
			// type S* suffix
			i ^= -1;
			if (i >= MASK_N) {
				i ^= MASK_N;
				tmp_int2 = 0;
				sorted_suffix[s] ^= MASK_N;
				// for the last two round if there is no recursion
			}
			if (last_sstar_deleg == equal_deleg) {
				--new_alphabet_size;
				new_pattern[--j] = ~i;
			}
			else {
				last_sstar_deleg = equal_deleg;
				new_pattern[--j] = i;
			}
			// temporarily use new_pattern to store the
			// list of s where type[sorted_suffix[s]] is 2
		}
		if (tmp_int2)
			equal_deleg = s - 1;
		else
			tmp_int2 = 1;
	}
	if (j != 0) {
		fprintf(stderr, "j = %d\n");
		assert(0);
	}

	int do_recursion = 1;

	// if there are no dup char in the new string then we're done.
	// if the new_alphabet_size <= 1 ($ and at most one different char)
	// this two check cover all new_length <= 2 case.
	if (new_alphabet_size == new_length || new_alphabet_size == 1)
		do_recursion = 0;

if (do_recursion) {
	gettimeofday(&tv, NULL);
	fprintf(stderr, " %d prepare rec at %d.%06d\n",
		length, tv.tv_sec, tv.tv_usec);
	// sorted_suffix is free to use. We store the sorted substring
	// value for each substring start at position k in
	// sorted_suffix[k / 2]
	/* generate the new pattern string and compute the
	 * new_bucket_end, new_bucket_p */
	int *new_bucket_end =
		(int *)_suffix_array_malloc(new_alphabet_size * sizeof(int));
	int *new_bucket_p =
		(int *)_suffix_array_malloc(new_alphabet_size * sizeof(int));
	j = 0;
	for (i = 0; i < new_length;) {
		s = new_pattern[i];
		++i;
		if (s < 0)
			sorted_suffix[(~s) >> 1] = j;
		else {
			sorted_suffix[s >> 1] = j;
			new_bucket_end[j] = i;
			++j;
		}
	}
	if (j != new_alphabet_size) {
		fprintf(stderr, "j = %d\n");
		assert(0);
	}
	gettimeofday(&tv, NULL);
	fprintf(stderr, " %d prepare rec. at %d.%06d\n",
		length, tv.tv_sec, tv.tv_usec);

	/* Generate the new pattern string and then find all type S*
	 * positions for the new pattern and put them into their
	 * bucket */
	gettimeofday(&tv, NULL);
	fprintf(stderr, "init %d at %d.%06d\n", new_length,
		tv.tv_sec, tv.tv_usec);

	int new_new_length = 0;
	int *new_new_orig_pos =
		(int *)_suffix_array_malloc(new_length * sizeof(int)); // FIXME
	memset(new_sorted_suffix, 0, new_length * sizeof(int));
	memcpy(new_bucket_p, new_bucket_end, new_alphabet_size * sizeof(int));
	int same = 0;
	int this_chr = -1;
	c0 = 0;
	b = new_sorted_suffix + new_bucket_p[c0];
	char r = 0;
	for (i = 0; i < new_length; ++i) {
		c1 = this_chr;
		this_chr = sorted_suffix[new_orig_pos[i] >> 1];
		new_pattern[i] = this_chr;
		if (c1 < this_chr) {
			if (r == 1) {
				new_new_orig_pos[new_new_length++] = same;
				if (c1 != c0) {
					new_bucket_p[c0] =
						b - new_sorted_suffix;
					b = new_sorted_suffix
						+ new_bucket_p[c0 = c1];
				}
				*(--b) = ~same ^ MASK_N;
				r = 0;
			}
			// same = i; // can be omitted
		}
		else if (c1 > this_chr) {
			r = 1;
			same = i;
		}
	}
	new_bucket_p[c0] = b - new_sorted_suffix;

	// new_alphabet_size is at least 2
	for (i = new_alphabet_size - 1; i >= 0; --i)
		if ((tmp_int1 = new_bucket_p[i]) < new_bucket_end[i])
			new_sorted_suffix[tmp_int1] ^= MASK_N;

	/* recursion */
	_suffix_array_worker<int>(new_pattern, new_length, new_new_length,
				  new_new_orig_pos, new_alphabet_size,
				  new_sorted_suffix,
				  new_bucket_end, new_bucket_p);

	gettimeofday(&tv, NULL);
	fprintf(stderr, " %d post rec start at %d.%06d\n",
		length, tv.tv_sec, tv.tv_usec);
	/* fix the position of the S* suffixes */
	c0 = alphabet_size;
	b = sorted_suffix + length;
	int *bb;
	for (i = new_length - 1; i >= 0; --i) {
		s = new_orig_pos[new_sorted_suffix[i]];
		c1 = pattern[s];
		if (c1 < c0) {
			bb = sorted_suffix + bucket_end[c0 = c1];
			while (b > bb)
				*(--b) = 0;
			if (c0 == 0)
				bb = sorted_suffix;
			else
				bb = sorted_suffix + bucket_end[c0 - 1];
		}
		*(--b) = ~s;
	}
	while (b > sorted_suffix)
		*(--b) = 0;

	gettimeofday(&tv, NULL);
	fprintf(stderr, " %d post rec finish at %d.%06d\n",
		length, tv.tv_sec, tv.tv_usec);
}

	/* sort */
	// first round
	// if |sorted_suffix[s]| - 1 is type S, save + in
	// sorted_suffix[s], otherwise save - in it.
	memcpy(bucket_p + 1, bucket_end, (alphabet_size - 1) * sizeof(int));
	i = length - 1;
	c0 = pattern[i];
	b = sorted_suffix + bucket_p[c0];
	if (pattern[i - 1] < c0)
		*b = i;
	else
		*b = ~i;
	++b;
	for (s = 0; s < length; ++s) {
		i = sorted_suffix[s];
		if (i < -1) {
			i = (~i) - 1;
			if ((c1 = pattern[i]) != c0) {
				bucket_p[c0] = b - sorted_suffix;
				b = sorted_suffix + bucket_p[c0 = c1];
			}
			if (i > 0 && pattern[i - 1] < c0)
				*b = i;
			else
				*b = ~i;
			++b;
		}
	}
	// second round
 	memcpy(bucket_p, bucket_end, alphabet_size * sizeof(int));
	b = sorted_suffix + bucket_p[c0];
	for (s = length - 1; s >= 0; --s) {
		i = sorted_suffix[s];
		if (i > 0) {
			--i;
			if ((c1 = pattern[i]) != c0) {
				bucket_p[c0] = b - sorted_suffix;
				b = sorted_suffix + bucket_p[c0 = c1];
			}
			--b;
			if (i > 0 && pattern[i - 1] <= c0)
				*b = i;
			else
				*b = ~i;
		}
		else
			sorted_suffix[s] ^= -1;
	}
	gettimeofday(&tv, NULL);
	fprintf(stderr, "end of %d at %d.%06d\n", length, tv.tv_sec, tv.tv_usec);
}

void suffix_array(char *pattern, int length) {
	int alphabet_size = 256;
	int *sorted_suffix = (int *)_suffix_array_malloc(length * sizeof(int));
	int *bucket_end = (int *)_suffix_array_malloc(alphabet_size * sizeof(int));
	int *bucket_p = (int *)_suffix_array_malloc(alphabet_size * sizeof(int));
	int *new_orig_pos;

	assert(length > 2);

	int new_length;
	new_length = _suffix_array_init(pattern, length, alphabet_size,
					sorted_suffix,
					bucket_end, bucket_p, &new_orig_pos);
	_suffix_array_worker<char>(pattern, length, new_length, new_orig_pos,
				   alphabet_size, sorted_suffix,
				   bucket_end, bucket_p);

	int i;
	for (i = 0; i < length; ++i)
		printf("%d\n", sorted_suffix[i]);
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
