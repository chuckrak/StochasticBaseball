from collections import defaultdict
import math
import random
def max_packages_on_shelf(client):
    shelf = set()
    max_packages = 0
    processed_packages = 1
    processed = defaultdict(bool)
    for i in range(len(client)):
        # If the package can be collected by the first client in line, remove it from the shelf
        desired_package = client[i]
        # print(f"Desired package is {desired_package}")
        # print(f"Shelf before processing is {shelf}")
        if desired_package in shelf:
            shelf.remove(desired_package)
            processed[desired_package] = True
        
        else:
            # print(f"Processed package is {processed_packages}")
            for i in range(processed_packages, desired_package):
                if not processed[i]:
                    shelf.add(i)
                    processed_packages += 1
            # print(f"Shelf after processing is {shelf}")
            processed[desired_package] = True
        max_packages = max(max_packages, len(shelf))

    return max_packages

def max_packages_fast(client):
    max_packages = -math.inf
    for i in range(len(client)):
        packages_on_shelf = client[i] - (i + 1)
        max_packages = max(max_packages, packages_on_shelf)
    return max_packages





# Test cases
print(max_packages_on_shelf([3, 2, 4, 5, 1]), max_packages_fast([3, 2, 4, 5, 1]))  # Output: 2
print(max_packages_on_shelf([1, 2, 3, 4, 5]), max_packages_fast([1, 2, 3, 4, 5]))  # Output: 0
print(max_packages_on_shelf([3, 2, 7, 5, 4, 1, 6]), max_packages_fast([3, 2, 7, 5, 4, 1, 6]))  # Output: 4
print(max_packages_on_shelf([5, 3, 7, 4, 6, 2, 1]), max_packages_fast([5, 3, 7, 4, 6, 2, 1]))  # Output: 4


rand_test_case = [i for i in range(1, 100000)]

for i in range(10):
    rand = random.shuffle(rand_test_case)
    print(max_packages_on_shelf(rand_test_case), max_packages_fast(rand_test_case))