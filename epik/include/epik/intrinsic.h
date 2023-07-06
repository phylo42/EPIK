#ifndef EPIK_INTRINSIC_H
#define EPIK_INTRINSIC_H

#ifdef EPIK_SSE
#include <vector>
#include <xmmintrin.h>
#include <i2l/phylo_kmer.h>

template <class T>
void update_vector(std::vector<i2l::phylo_kmer::score_type>& vec,
                   std::vector<size_t>& counts,
                   std::vector<i2l::phylo_kmer::branch_type>& edges,
                   const T& updates) {
    int i = 0;

    // Process updates in blocks of 4 as long as possible
    for (; i <= (int)updates.size() - 4; i += 4) {
        // Load score updates
        __m128 simdValues = _mm_set_ps(updates[i].score, updates[i+1].score,
                                        updates[i+2].score, updates[i+3].score);

        // Load the current scores
        __m128 currentValues = _mm_set_ps(vec[updates[i].branch], vec[updates[i+1].branch],
                                          vec[updates[i+2].branch], vec[updates[i+3].branch]);

        // SIMD Add
        __m128 newValues = _mm_add_ps(currentValues, simdValues);

        // Store the new values back in the vector individually
        for (int j = 0; j < 4; j++) {
            vec[updates[i+j].branch] = newValues[j];
            if (counts[updates[i+j].branch] == 0)
            {
                edges.push_back(updates[i+j].branch);
            }
            counts[updates[i+j].branch]++;
        }
    }

    // Process the remaining updates
    for (; i < (int)updates.size(); i++) {
        vec[updates[i].branch] += updates[i].score;
        counts[updates[i].branch]++;
    }
}
#endif

//#define EPIK_AVX2
#ifdef EPIK_AVX2

#include <vector>
#include <immintrin.h>
#include <i2l/phylo_kmer.h>

template <class T>
void update_vector(std::vector<i2l::phylo_kmer::score_type>& vec,
                   std::vector<size_t>& counts,
                   std::vector<i2l::phylo_kmer::branch_type>& edges,
                   const T& updates) {

    int i = 0;

    // 256 bits / 32-bit float = 8 floats
    constexpr int simdWidth = 8;
    int indices[simdWidth];

    // Process updates in blocks of simdWidth as long as possible
    for (; i <= (int)updates.size() - simdWidth; i += simdWidth) {
        // Prepare indices
        for (int j = 0; j < simdWidth; j++) {
            indices[j] = updates[i + j].branch;
        }

        // Load score updates
        __m256 simdValues = _mm256_set_ps(updates[i].score, updates[i+1].score, updates[i+2].score,
                                          updates[i+3].score, updates[i+4].score, updates[i+5].score,
                                          updates[i+6].score, updates[i+7].score);

        // Load indices and gather the current scores
        __m256i simdIndices = _mm256_loadu_si256((const __m256i*)indices);
        __m256 currentValues = _mm256_i32gather_ps(vec.data(), simdIndices, 4);

        // SIMD Add
        __m256 newValues = _mm256_add_ps(currentValues, simdValues);

        // Store the new values back in the vector individually
        float results[simdWidth];
        _mm256_storeu_ps(results, newValues);

        for (int j = 0; j < simdWidth; j++) {
            vec[updates[i+j].branch] = results[j];
            if (counts[updates[i+j].branch] == 0)
            {
                edges.push_back(updates[i+j].branch);
            }
            counts[updates[i+j].branch]++;
        }
    }

    // Process the remaining updates
    for (; i < (int)updates.size(); i++) {
        vec[updates[i].branch] += updates[i].score;
        if (counts[updates[i].branch] == 0)
        {
            edges.push_back(updates[i].branch);
        }
        counts[updates[i].branch]++;
    }
}

#endif


#ifdef EPIK_AVX512
#include <vector>
#include <immintrin.h>
#include <i2l/phylo_kmer.h>

template <class T>
void update_vector(std::vector<i2l::phylo_kmer::score_type>& vec,
                   std::vector<size_t>& counts,
                   std::vector<i2l::phylo_kmer::branch_type>& edges,
                   const T& updates) {

    int i = 0;
    constexpr int simdWidth = 16;  // 512 bits / 32-bit float = 16 floats
    int indices[simdWidth];

    // Process updates in blocks of simdWidth as long as possible
    for (; i <= (int)updates.size() - simdWidth; i += simdWidth) {
        // Prepare indices
        for (int j = 0; j < simdWidth; j++) {
            indices[j] = updates[i + j].branch;
        }

        // Load score updates
        __m512 simdValues = _mm512_set_ps(updates[i].score, updates[i+1].score, updates[i+2].score,
                                          updates[i+3].score, updates[i+4].score, updates[i+5].score,
                                          updates[i+6].score, updates[i+7].score, updates[i+8].score,
                                          updates[i+9].score, updates[i+10].score, updates[i+11].score,
                                          updates[i+12].score, updates[i+13].score, updates[i+14].score,
                                          updates[i+15].score);

        // Load indices and gather the current scores
        __m512i simdIndices = _mm512_loadu_si512((__m512i const*)indices);
        __m512 currentValues = _mm512_i32gather_ps(simdIndices, vec.data(), 4);

        // SIMD Add
        __m512 newValues = _mm512_add_ps(currentValues, simdValues);

        // Store the new values back in the vector individually
        float results[simdWidth];
        _mm512_storeu_ps(results, newValues);

        for (int j = 0; j < simdWidth; j++) {
            vec[updates[i+j].branch] = results[j];
            if (counts[updates[i+j].branch] == 0)
            {
                edges.push_back(updates[i+j].branch);
            }
            counts[updates[i+j].branch]++;
        }
    }

    // Process the remaining updates
    for (; i < (int)updates.size(); i++) {
        vec[updates[i].branch] += updates[i].score;
        if (counts[updates[i].branch] == 0)
        {
            edges.push_back(updates[i].branch);
        }
        counts[updates[i].branch]++;
    }
}


#endif

#endif
