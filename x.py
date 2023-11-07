###There is an array A of N non-negative integers. Any two initial elements of A that are adjacent can be replaced with their merged equivalent. 
# For example, given A = [2, 3, 15], pair (2, 3) can be replaced with 23, resulting in array (23, 15], 
# and pair (3, 15) can be replaced with 315, resulting in array [2, 315].
# The result of the merge cannot be merged any further, so we can't get 2315 in the example above. 
# What is the maximum possible sum of elements of A after any number of merges?
# Write a function:
# def solution (A)
# that, given an array A of N non-negative integers, returns the maximum sum of
# elements of A after any number of merges.
# def merge_elements(a, b):
#     return int(str(a) + str(b))

# def dp_merge_sum(A, dp, i):
#     # Base case: If the current index is out of bounds, return 0
#     if i >= len(A):
#         return 0

#     # If the maximum sum for this index has already been computed, return it
#     if dp[i] is not None:
#         return dp[i]

#     # Consider the three possible options: merge with next element, merge with the element before, or do not merge
#     merge_with_next = A[i]
#     merge_with_prev = A[i]
#     if i + 1 < len(A):
#         merge_with_next = merge_elements(A[i], A[i + 1]) + dp_merge_sum(A, dp, i + 2)
#     if i - 1 >= 0:
#         merge_with_prev = merge_elements(A[i - 1], A[i]) + dp_merge_sum(A, dp, i + 1)
    
#     not_merged = A[i] + dp_merge_sum(A, dp, i + 1)

#     # Choose the maximum sum for the current position and store it in the dp table
#     dp[i] = max(merge_with_next, merge_with_prev, not_merged)
#     return dp[i]

# def solution(A):

#     # Initialize the dp table with None (indicating that values are not computed yet)
#     dp = [None] * len(A)

#     # Find the maximum possible sum starting from the first element (index 0)
#     max_sum = dp_merge_sum(A, dp, 0)

#     return max_sum

# # Test cases
# print(solution([2, 2, 3, 5, 4, 0]))  # Output: 97
# print(solution([3, 19, 191, 91, 3]))  # Output: 20107
# print(solution([12, 6, 18, 10, 1, 0]))  # Output: 1946
# print(solution([2, 1, 0, 1, 2, 9, 1, 0]))  # Output: 124


# Examples:
# 1. Given A = [2, 2, 3, 5, 4, 0] the function should return 97. We can merge elements of the following pairs: (2, 2), (3, 5) and (4, 0). This results in A = [22, 35, 40], which
# sums up to 97.
# 2. Given A = [3, 19, 191, 91, 3] the function should return 20107. We can merge elements of the following pairs: (19, 191) and (91, 3). This results in A = [3,
# 19191, 913], which sums up to 20107.
# 3. Given A = [12, 6, 18, 10, 1, 0], the function should return 1946. The merges should make A = [126, 1810, 10], which sums up to 1946.
# 4. Given A = [2, 1, 0, 1, 2, 9, 1, O], the function should return 124. The merges should make A = [21, 0, 12, 91, 0], which sums up to 124




def merge_elements(a, b):
    return int(str(a) + str(b))

def solution(A):
    n = len(A)

    # Initialize the dynamic programming table
    dp = [0] * (n + 1)

    # Fill in the dynamic programming table from right to left
    for i in range(n - 1, -1, -1):
        # Option 1: Merge with the next element (if possible)
        merged_sum = 0
        if i + 1 < n:
            merged_sum = merge_elements(A[i], A[i + 1]) + dp[i + 2]

        # Option 2: Do not merge anything
        not_merged_sum = A[i] + dp[i + 1]

        # Choose the maximum sum for the current index
        dp[i] = max(merged_sum, not_merged_sum)

    # The maximum possible sum will be stored in dp[0]
    return dp[0]

# Test cases
print(solution([2, 3, 15]))
print(solution([2, 2, 3, 5, 4, 0]))  # Output: 97
print(solution([3, 19, 191, 91, 3]))  # Output: 20107
print(solution([12, 6, 18, 10, 1, 0]))  # Output: 1946
print(solution([2, 1, 0, 1, 2, 9, 1, 0]))  # Output: 124
