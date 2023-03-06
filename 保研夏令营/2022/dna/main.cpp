#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <stack>
#include <algorithm>

#define CORE_POINT 1
#define BORDER_POINT 2
#define OUTLIER_POINT 3

/* 读取文件 */
void read_file(const std::string& path, std::string &file_str){
    std::ifstream in(path);
    if(in.fail())exit(0);
    std::ostringstream stream;
    stream << in.rdbuf();
    file_str = stream.str();
    in.close();
}

/* 三个数中最小的 */
uint32_t min_of_3(uint32_t val1, uint32_t val2, uint32_t val3){
    std::vector<uint32_t> vec = {val1, val2, val3};
    uint32_t min = UINT32_MAX;
    for(auto val : vec)
        if(val < min)
            min = val;
    return min;
}

void LoadData(std::string &data_file_path, uint32_t &n, uint32_t &e, uint32_t&p, std::vector<std::string> &dna_seqs){
    std::string file_content;
    read_file(data_file_path, file_content);
    std::stringstream content_stream(file_content);
    content_stream >> n >> e >> p;
    std::string line;
    content_stream >> line;
    while (!content_stream.eof()){
        if(!line.empty())
            dna_seqs.push_back(line);
        content_stream >> line;
    }
}

uint32_t LevDist(const std::string &a, const std::string &b){
    uint32_t a_length = a.size(), b_length = b.size();
    std::vector<std::vector<uint32_t>> dist_vec(a_length + 1, std::vector<uint32_t>(b_length + 1, 0));
    for(int i = 0; i <= a_length; i++)
        dist_vec[i][0] = i;
    for(int j = 0; j <= b_length; j++)
        dist_vec[0][j] = j;
    for(int i = 1; i <= a_length; i++)
        for(int j = 1; j <= b_length; j++){
            uint32_t val_1 = dist_vec[i-1][j]+1;
            uint32_t val_2 = dist_vec[i][j-1]+1;
            uint32_t val_3 = dist_vec[i-1][j-1];
            if(a[i-1] != b[j-1])
                val_3+=1;
            dist_vec[i][j] = min_of_3(val_1, val_2, val_3);
        }
    return dist_vec[a_length][b_length];
}

std::vector<std::string> RangeQuery(const std::vector<std::string> &dna_set, const std::string &target_dna, uint32_t eps){
    std::vector<std::string> less_than_eps_dna_set;
    for(const auto& dna_seq:dna_set){
        uint32_t dist = LevDist(dna_seq, target_dna);
        if(dist <= eps)
            less_than_eps_dna_set.push_back(dna_seq);
    }
    return less_than_eps_dna_set;
}

std::vector<std::vector<std::string>> DBSCAN(
        const std::vector<std::string> &dna_set,
        uint32_t eps, uint32_t minpts, uint32_t &num_cores, uint32_t &num_borders,
        uint32_t &num_outliers, uint32_t &num_clusters,
        std::vector<uint32_t> &point_class, std::vector<std::vector<uint32_t>> &dist_vec){
    uint32_t dna_num = dna_set.size();
    // 记录每个点在其半径为 eps 的圆形范围内的点
    std::vector<std::vector<uint32_t>> less_than_eps_vec = std::vector<std::vector<uint32_t>>(dna_num, std::vector<uint32_t>());
    // 记录每个点在图上的临接点
    std::vector<std::set<uint32_t>> near_point_set = std::vector<std::set<uint32_t>>(dna_num, std::set<uint32_t>());
    std::vector<bool> point_reached = std::vector<bool>(dna_num, false);
    std::vector<std::vector<std::string>> clusters;

    // STEP1: 计算距离, 并记录在其半径为 eps 的圆形范围内的点
    for(uint32_t i = 0; i < dna_num; i++)
        for(uint32_t j = 0; j < dna_num ; j++) {
            if(i == j)
                continue;
            uint32_t dist = LevDist(dna_set[i], dna_set[j]);
            dist_vec[i][j] = dist;
            if(dist <= eps)
                less_than_eps_vec[i].push_back(j);
        }

    // STEP2: 判断点的类型
    for(uint32_t i = 0; i < dna_num; i++)
        if(less_than_eps_vec[i].size() + 1 >= minpts){
            point_class[i] = CORE_POINT;
            for(auto point : less_than_eps_vec[i]) {
                near_point_set[point].insert(i);
                near_point_set[i].insert(point);
                if (point_class[point] != CORE_POINT)
                    point_class[point] = BORDER_POINT;
            }
        }

    // STEP3: 计数
    for(uint32_t i = 0; i < dna_num; i++){
        if(point_class[i] == CORE_POINT)
            num_cores ++;
        else if(point_class[i] == BORDER_POINT)
            num_borders ++;
        else
            num_outliers ++;
    }

    // STEP4: 使用深度优先搜索进行聚类
    for(int i = 0; i < dna_num; i++)
        if(point_class[i] != OUTLIER_POINT && !point_reached[i]){
            std::vector<std::string> cluster;
            std::stack<uint32_t> point_stack;
            point_stack.push(i);
            point_reached[i] = true;
            cluster.push_back(dna_set[i]);
            while(!point_stack.empty()){
                uint32_t point = point_stack.top();
                point_stack.pop();
                for(auto near_point : near_point_set[point])
                    if(!point_reached[near_point]){
                        point_reached[near_point] = true;
                        point_stack.push(near_point);
                        cluster.push_back(dna_set[near_point]);
                    }
            }
            clusters.push_back(cluster);
        }
    num_clusters = clusters.size();
    return clusters;
}

uint32_t MinEPS(const std::vector<std::string> &dna_set, uint32_t minpts,
                uint32_t prev_eps,const std::vector<uint32_t> &point_class,
                const std::vector<std::vector<uint32_t>> &dist_vec){
    //记录擦除掉每个OUTLIER_POINT需要的最小eps
    std::vector<uint32_t> new_eps_vec;
    for(int i = 0; i < dna_set.size(); i++)
        if(point_class[i] == OUTLIER_POINT){
            uint32_t new_min_eps_of_point;
            //该OUTLIER_POINT作为CORE_POINT的情况
            auto dist_vec_from_i = dist_vec[i];
            std::sort(dist_vec_from_i.begin(), dist_vec_from_i.end());
            uint32_t new_min_eps_core_i = dist_vec_from_i[minpts - 1];
            new_min_eps_of_point = new_min_eps_core_i;

            //该OUTLIER_POINT作为BORDER_POINT的情况
            for(int j = 0; j < dna_set.size(); j++){
                auto dist_vec_from_j = dist_vec[j];
                std::sort(dist_vec_from_j.begin(), dist_vec_from_j.end());
                uint32_t new_min_eps_core_j = std::max(dist_vec_from_j[minpts - 1], dist_vec[j][i]);
                if(new_min_eps_of_point > new_min_eps_core_j) new_min_eps_of_point = new_min_eps_core_j;
            }

            //现在已经找到该点的最小eps
            new_eps_vec.push_back(new_min_eps_of_point);
        }
    //返回全部OUTLIER_POINT的最小eps中最大的那个
    std::sort(new_eps_vec.begin(), new_eps_vec.end());
    //针对无OUTLIER_POINT的情况
    if(new_eps_vec.empty()) return prev_eps;
    return new_eps_vec[new_eps_vec.size()-1];
}

int main() {
    std::string data_file_path;
    uint32_t n,e,p, num_cores = 0, num_borders = 0, num_outliers = 0, num_clusters = 0;
    std::vector<std::string> dna_seqs;
    std::vector<uint32_t> point_class;
    std::vector<std::vector<uint32_t>> dist_vec;

    /* 2.1 */
    std::cin >> data_file_path;
    LoadData(data_file_path, n, e, p, dna_seqs);
    uint32_t min_length = UINT32_MAX, max_length = 0;
    for(const auto& dna_seq : dna_seqs){
        uint32_t seq_size = dna_seq.size();
        if(seq_size > max_length)
            max_length = seq_size;
        if(seq_size < min_length)
            min_length = seq_size;
    }
    std::cout <<min_length << " " << max_length <<std::endl;

    /* 2.2 */
    if(dna_seqs.size()<2)exit(0);
    std::cout << LevDist(dna_seqs[0], dna_seqs[1]) << std::endl;

    /* 2.3 */
    std::cout << RangeQuery(dna_seqs, dna_seqs[0], e).size() <<std::endl;

    /* 2.4 */
    uint32_t dna_num = dna_seqs.size();
    point_class = std::vector<uint32_t>(dna_num, OUTLIER_POINT);
    dist_vec = std::vector<std::vector<uint32_t>>(dna_num, std::vector<uint32_t>(dna_num, 0));
    auto clusters = DBSCAN(dna_seqs, e, p, num_cores, num_borders, num_outliers, num_clusters, point_class, dist_vec);
    std::cout <<num_cores << " "<< num_borders << " " << num_outliers <<" "<< num_clusters << " "<<std::endl;

    /* 2.5 */
    std::cout << MinEPS(dna_seqs, p, e, point_class, dist_vec) << std::endl;
}
