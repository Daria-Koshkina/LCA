#include <iostream>
#include <vector>

const int INDEXSHIFT = 1;

class Lowest_common_ancestor {
 private:
    int vertices_count;
    std::vector<std::vector<int>> tree;

    int block_size;
    int block_count;
    std::vector<int> first;
    std::vector<int> order;
    std::vector<int> depths;
    std::vector<int> logarithms_2;
    std::vector<std::vector<int>> sparse_table;
    std::vector<std::vector<std::vector<int>>> blocks;
    std::vector<int> block_mask;

public:
    explicit Lowest_common_ancestor(const std::vector<int>& vertices) {
        vertices_count = vertices.size();
        build_tree(vertices);
        preprocess();
    }

    // gets lca for two vertex
    int get_lowest_common_ancestor(const int first_vertex, const int second_vertex);

private:
    // does dfs in which additionally remembers the order
    // of traversing the vertices and their heights
    void depth_first_search(
        const int vertex,
        const int parent,
        const int height);

    // gets index of vertex with lowest height
    int get_minimum_by_height(
        const int first_vertex,
        const int second_vertex);

    // builds tree of given vertices like adjency list
    void build_tree(const std::vector<int>& vertices);

    // does all preprocessing of lca and rmq
    void preprocess();

    // get lca in one block
    int lowest_common_ancestor_in_block(
        const int block,
        const int left_vertex,
        const int right_vertex);
};

// containes input data
struct Input_data {
    std::vector<int> vertices;
    int query_count;
    int64_t x_coefficient;
    int64_t y_coefficient;
    int64_t z_coefficient;
    int first_query_vertex;
    int second_query_vertex;
};

// reads input data
Input_data get_input_data(std::istream& input);

// gets answer on queries
int64_t get_sum_of_query_answers(const Input_data& inputDataSet);

// prints result
void print_result(const int64_t sum, std::ostream& output);


void Lowest_common_ancestor::depth_first_search(
    const int vertex,
    const int parent,
    const int height) {
    first[vertex] = order.size();
    order.push_back(vertex);
    depths[vertex] = height;

    for (int next_vertex : tree[vertex]) {
        if (next_vertex == parent) {
            continue;
        }
        depth_first_search(next_vertex, vertex, height + 1);
        order.push_back(vertex);
    }
}

int Lowest_common_ancestor::get_minimum_by_height(
    const int first_vertex,
    const int second_vertex) {
    if (depths[order[first_vertex]] < depths[order[second_vertex]]) {
        return first_vertex;
    }
    else {
        return second_vertex;
    }
}

void Lowest_common_ancestor::build_tree(const std::vector<int>& vertices) {
    tree.resize(vertices.size());
    for (int i = 1; i <= vertices.size() - 1; ++i) {
        tree[vertices[i - INDEXSHIFT]].push_back(i);
        tree[i].push_back(vertices[i - INDEXSHIFT]);
    }
}

void Lowest_common_ancestor::preprocess() {
    first.assign(vertices_count, -1);
    depths.assign(vertices_count, 0);
    order.reserve(2 * vertices_count);
    depth_first_search(0, -1, 0);

    // precompute all log values
    int order_vertices_count = order.size();
    logarithms_2.resize(0);
    logarithms_2.push_back(-1);
    for (int i = 1; i <= order_vertices_count; i++) {
        logarithms_2.push_back(logarithms_2[i / 2] + 1);
    }

    block_size = std::max(1, logarithms_2[order_vertices_count] / 2);
    block_count = (order_vertices_count + block_size - 1) / block_size;

    // build sparse table
    sparse_table.assign(block_count,
        std::vector<int>(logarithms_2[block_count] + 1));
    for (int i = 0, j = 0, b = 0; i < order_vertices_count; i++, j++) {
        if (j == block_size) {
            j = 0, b++;
        }
        if (j == 0 || get_minimum_by_height(i, sparse_table[b][0]) == i) {
            sparse_table[b][0] = i;
        }
    }
    for (int l = 1; l <= logarithms_2[block_count]; l++) {
        for (int i = 0; i < block_count; i++) {
            int next_index = i + (1 << (l - 1));
            if (next_index >= block_count) {
                sparse_table[i][l] = sparse_table[i][l - 1];
            }
            else {
                sparse_table[i][l] =
                    get_minimum_by_height(
                        sparse_table[i][l - 1],
                        sparse_table[next_index][l - 1]);
            }
        }
    }

    // get mask for each block
    block_mask.assign(block_count, 0);
    for (int i = 0, j = 0, b = 0; i < order_vertices_count; i++, j++) {
        if (j == block_size)
            j = 0, b++;
        if (j > 0 && (i >= order_vertices_count ||
            get_minimum_by_height(i - 1, i) == i - 1))
            block_mask[b] += 1 << (j - 1);
    }

    // preprocess RMQ for each unique block
    int rmq_blocks_count = 1 << (block_size - 1);
    blocks.resize(rmq_blocks_count);
    for (int block_index = 0; block_index < block_count; block_index++) {
        int mask = block_mask[block_index];
        if (!blocks[mask].empty()) {
            continue;
        }
        blocks[mask].assign(block_size, std::vector<int>(block_size));
        for (int left = 0; left < block_size; left++) {
            blocks[mask][left][left] = left;
            for (int right = left + 1; right < block_size; right++) {
                blocks[mask][left][right] = blocks[mask][left][right - 1];
                if (block_index * block_size + right < order_vertices_count) {
                    blocks[mask][left][right] =
                        get_minimum_by_height(
                            block_index * block_size +
                            blocks[mask][left][right],
                            block_index * block_size + right)
                        - block_index * block_size;
                }
            }
        }
    }
}

int Lowest_common_ancestor::lowest_common_ancestor_in_block(
    const int block,
    const int left_vertex,
    const int right_vertex) {
    return blocks[block_mask[block]][left_vertex][right_vertex] +
        block * block_size;
}

int Lowest_common_ancestor::get_lowest_common_ancestor(
    const int first_vertex,
    const int second_vertex) {
    int left_vertex = first[first_vertex];
    int right_vertex = first[second_vertex];
    if (left_vertex > right_vertex) {
        std::swap(left_vertex, right_vertex);
    }
    int left_block = left_vertex / block_size;
    int right_block = right_vertex / block_size;
    if (left_block == right_block) {
        int64_t answer = lowest_common_ancestor_in_block(
            left_block,
            left_vertex % block_size,
            right_vertex % block_size);
        return order[answer];
    }
    int answer1 = lowest_common_ancestor_in_block(
        left_block,
        left_vertex % block_size,
        block_size - 1);
    int answer2 = lowest_common_ancestor_in_block(
        right_block,
        0,
        right_vertex % block_size);
    int answer = get_minimum_by_height(answer1, answer2);
    if (left_block + 1 < right_block) {
        int logarithm = logarithms_2[right_block - left_block - 1];
        int answer3 = sparse_table[left_block + 1][logarithm];
        int answer4 = sparse_table[right_block - (1 << logarithm)][logarithm];
        answer = get_minimum_by_height(
            answer,
            get_minimum_by_height(answer3, answer4));
    }
    return order[answer];
}

Input_data get_input_data(std::istream& input = std::cin) {
    Input_data inputDataSet;
    int vertexCount, query_count;
    input >> vertexCount >> query_count;
    inputDataSet.query_count = query_count;
    inputDataSet.vertices = std::vector<int>(vertexCount);
    for (int i = 0; i < vertexCount - INDEXSHIFT; ++i) {
        input >> inputDataSet.vertices[i];
    }
    input >> inputDataSet.first_query_vertex >> inputDataSet.second_query_vertex;
    input >> inputDataSet.x_coefficient >> inputDataSet.y_coefficient >> inputDataSet.z_coefficient;
    return inputDataSet;
}

int64_t get_sum_of_query_answers(const Input_data& inputDataSet) {
    Lowest_common_ancestor lcaTree = Lowest_common_ancestor(inputDataSet.vertices);
    int64_t sum = 0;
    int lcaResult = 0;
    int first_query_vertex = inputDataSet.first_query_vertex;
    int second_query_vertex = inputDataSet.second_query_vertex;
    const int vertexCount = inputDataSet.vertices.size();
    const int64_t x_coefficient = inputDataSet.x_coefficient;
    const int64_t y_coefficient = inputDataSet.y_coefficient;
    const int64_t z_coefficient = inputDataSet.z_coefficient;
    for (int i = 0; i < inputDataSet.query_count; ++i) {
        int oneVertex = (first_query_vertex + lcaResult) % vertexCount;
        lcaResult = lcaTree.get_lowest_common_ancestor(oneVertex, second_query_vertex);
        sum += lcaResult;
        first_query_vertex = (x_coefficient * first_query_vertex + y_coefficient * second_query_vertex + z_coefficient)
            % vertexCount;
        second_query_vertex = (x_coefficient * second_query_vertex + y_coefficient * first_query_vertex + z_coefficient)
            % vertexCount;
    }
    return sum;
}

void print_result(const int64_t sum, std::ostream& output = std::cout) {
    output << sum;
}

int main()
{
    Input_data startData = get_input_data();
    const int64_t result = get_sum_of_query_answers(startData);
    print_result(result);

    return 0;
}
