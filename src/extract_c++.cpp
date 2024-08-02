//version 2.0


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <nlohmann/json.hpp>
#include <Eigen/Dense>
#include <limits>
#include <algorithm>
// 使用命名空间简化代码
using json = nlohmann::json;
using Eigen::MatrixXd;
using Eigen::Vector2d;



// 计算两个点之间的距离
double distance(const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
    return (a - b).norm();
}

//Graham扫描算法求凸包，基本原理是利用叉积判断两个向量的走向是否是顺时针or逆时针（前一个向量到后一个向量是右拐还是左拐）
        //叉积是如何判断左右拐的？根据右手定则，若向量右拐，则叉积方向沿z轴负向，反之正向，而在二维平面中可以视作叉乘结果小于0右拐，大于0左拐
    /*Graham扫描执行流程：
        （1）先把所有点按照横坐标进行排序，从小到大
        （2）随后，从最小的点开始，按照横坐标次序依次放入栈中，并随时判断栈顶的三个点是否依旧构成连续右转的关系，如果出现左转的情况，设此时新的要加入点是C，
        C还没有入栈，但是判断出目前栈顶的两个点和C出现了左转关系，则开始出栈，直到栈顶的两个点和点C构成连续右转关系，停止出栈，将C入栈
        （3）重复第二步，从左到右找到上凸包，再从右到左找到下凸包，然后合成
                                    
    */
// 计算二维向量的叉积
double cross(const Vector2d& a, const Vector2d& b) {
    return a.x() * b.y() - a.y() * b.x();
}

//分别求上下凸包并合成
std::vector<Eigen::Vector2d> compute_convex_hull(std::vector<Eigen::Vector2d> points) { //利用Eigen库中用于表示二维向量的vector2d类

// 按 x 坐标排序，如果 x 坐标相同则按 y 坐标排序
    std::sort(points.begin(), points.end(), [](const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) { //使用标准库sort函数
        return (p1.x() < p2.x()) || (p1.x() == p2.x() && p1.y() < p2.y()); //排序准则是按x从小到大排序，如果 x 坐标相同则按 y 坐标排序
    });

    // 构建凸包的下半部分
    std::vector<Eigen::Vector2d> hull;//二维点栈hull
    for (const auto& point : points) {//遍历二维点集points里面的每一个点point
        while (hull.size() >= 2 && cross(hull[hull.size() - 1] - hull[hull.size() - 2], point - hull[hull.size() - 1]) <= 0) {
            //当栈hull的元素个数大于2个，但出现三个点的向量出现右转时出栈
            hull.pop_back(); 
        }
        hull.push_back(point); //一般情况一直入栈
    }

    size_t lower_hull_size = hull.size();//保存下半部分大小

    // 构建凸包的上半部分
    
    for (auto it = points.rbegin(); it != points.rend(); ++it) { //使用反向迭代器，从右往左遍历point
        while (hull.size() > lower_hull_size && cross(hull[hull.size() - 1] - hull[hull.size() - 2], *it - hull[hull.size() - 1]) <= 0) {
            //最右侧的起始点已经在栈内，而在上半部分加入第一个点时不满足hull的大小大于lower_hull_size，因此只有在加入上半部分第二个点时才会触发条件hull.size() > lower_hull_size
            //这个时候算上起始点，第一个点，已经满足了等效hull_size_higher>=2，所以hull.size() > lower_hull_size等效于higher_hull_size>=2,
            //后面使用迭代器指针代表点point
            hull.pop_back();
        }
        hull.push_back(*it);
    }

    hull.pop_back(); // 删除最后一个点，因为它与第一个点相同
    return hull;


}

//旋转卡壳法，不是点对应最远边，而是边对应最远点
//初步思路：先找主轴（长径），再找次主轴（短径）；其中轴不一定非得是点对点，也可以是点对边
//正因为这一点，所以我们在执行时，先找到距离最远的两点作为长径，然后在垂直长径的方向上从最左或最右逐点扫描到另一端距离（这个端可以是点也可以是两点间的边），找到最大距离作为短径
// 计算凸包的长径（最大距离）
std::pair<double, std::pair<Vector2d, Vector2d>> cal_major_axis(const std::vector<Eigen::Vector2d>& hull) {
    if (hull.size() <= 1) {
        std::cerr << "Convex hull has too few points!" << std::endl;
        return {0, {Vector2d(0, 0), Vector2d(0, 0)}};
    }
    if (hull.size() == 2) {
        return {distance(hull[0], hull[1]), {hull[0], hull[1]}};
    }

    double max_dist = 0;
    Vector2d p1, p2;
    for (size_t i = 0; i < hull.size(); ++i) {
        for (size_t j = i + 1; j < hull.size(); ++j) {
            double dis = distance(hull[i], hull[j]);
            if (dis > max_dist) {
                max_dist = dis;
                p1 = hull[i];
                p2 = hull[j];
            }
        }
    }
    return {max_dist, {p1, p2}};
}

// 计算点到直线的距离
double point_to_line_distance(const Vector2d& p, const Vector2d& a, const Vector2d& b) {
    return std::abs(cross(b - a, p - a) / (b - a).norm());
}
// 叉积是三角形面积，面积除以底边ab的模（norm()方法算模长）得到高，也就是点p到直线ab的距离
// abs（）负责求模



/*---初步排查是交点计算错误导致出现过大值---*/
//   ↑已解决 原因是垂线没有正常延伸，因而没有正确计算交点


//下个问题是为什么长径会不匹配，我认为长径应该至少是百分百匹配的
//短径不匹配的原因尚不清晰，有长有短



// 计算直线与线段的交点
Vector2d line_segment_intersection(const Vector2d& a1, const Vector2d& a2, const Vector2d& b1, const Vector2d& b2, bool& intersects) {
    Vector2d r = a2 - a1; //代表线段a1a2的方向
    Vector2d s = b2 - b1; //代表线段b1b2的方向
    double rxs = cross(r, s);
    
    if (rxs == 0) {
        // 线段平行或共线
        intersects = false;
        return Vector2d(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());
    }

    //有了方向，我们可以设两条线段的表达式
    //P(t)=a1+t*r=a1+t*(a2-a1)  线段a1a2
    //Q(u)=b1+u*s=b1+u*(b2-b1)  线段b1b2
    

    //然后我们可以求交点，将两个式子联立
    /*
    可以得到a1-b1=u*s-t*r，分解为两个标量方程
    (a1−b1)⋅s ⊥ =(u⋅s−t⋅r)⋅s ⊥    s⊥表示向量s的垂直向量（-sy,sx)
    因为r*s ⊥=rx*sy-ry*sx也就是r x s，这对于其他变量也是同理

    方程可化为cross(a1-b1,s)=u*cross(s,s)-t*cross(r,s)
    可以化简为：cross(a1-b1,s)=-t*cross(r,s)

    可以得到最终表达式↓
    */

    double t = cross(b1 - a1, s) / rxs;
    double u = cross(b1 - a1, r) / rxs;

    //通过这两个参数，我们可以确认两个线段是否相交
    //如果都在[0,1]之内，那么说明1倍线段长度，也就是a1a2和b1b2本身就相交
    //否则就不存在交点

    if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
        intersects = true;
        return a1 + t * r; //如果存在交点，就返回线段表达式
    } else {
        intersects = false;
        return Vector2d(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());
    }
}


// 计算凸包的短径（垂直于长径的最大长度）
double cal_minor_axis(const std::vector<Eigen::Vector2d>& hull, const Vector2d& p1, const Vector2d& p2) {
    Vector2d major_axis = p2 - p1;
    Vector2d normal(-major_axis.y(), major_axis.x()); //创建法向量normal
    normal.normalize(); // 归一化

    double max_length = 0;
    for (const auto& point : hull) {
        if (cross(major_axis, point - p1) > 0) {//叉积为正表示point-p1在major_axis左侧，用于筛选主轴一侧的点
            double max_distance = 0;
            for (size_t i = 0; i < hull.size(); ++i) {//遍历凸包所有边
                size_t next_i = (i + 1) % hull.size();
                if (cross(major_axis, hull[next_i] - p1) <= 0) {//判断目标边是否在右侧
                    bool intersects;
                    Vector2d intersection = line_segment_intersection(point - 1000 * normal, point + 1000 * normal, hull[i], hull[next_i], intersects);//计算交点，延长法向量长度确保相交
                    if (intersects) {
                        double dist = distance(point, intersection);
                        max_distance = std::max(max_distance, dist);
                    }
                }
            }
            max_length = std::max(max_length, max_distance);
        }
    }

    return max_length;
}



// 计算长径和短径的函数
std::pair<double, double> calculate_diameters(const std::vector<std::vector<double>>& points) {
    std::vector<Eigen::Vector2d> eigen_points;
    for (const auto& point : points) {
        eigen_points.emplace_back(point[0], point[1]);
    }
    std::vector<Eigen::Vector2d> hull = compute_convex_hull(eigen_points);
    auto [major_axis_length, major_axis_points] = cal_major_axis(hull);
    double minor_axis_length = cal_minor_axis(hull, major_axis_points.first, major_axis_points.second);
    return {major_axis_length, minor_axis_length};
}

int main() {
    // 打开 JSON 文件
    std::string file_path = "C:/code/extract_git/src/predict_ct_chest_vr-0722.json"; // 确保路径正确
    std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        return 1;
    }

    // 输出文件路径
    std::cout << "Reading JSON file from: " << file_path << std::endl;

    // 重新打开文件进行解析
    json data;
    try {
        file >> data;
    } catch (json::parse_error& e) {
        std::cerr << "JSON parse error: " << e.what() << std::endl;
        return 1;
    }

    // 提取 ct_nodule 数据
    if (!data.contains("ct_nodule")) {
        std::cerr << "JSON does not contain 'ct_nodule' key!" << std::endl;
        return 1;
    }

    auto ct_nodules = data["ct_nodule"];

    for (const auto& nodule : ct_nodules) {
        for (const auto& contour : nodule["contour3D"]) {
            std::vector<std::vector<double>> points;

            // 提取二维坐标点
            for (const auto& point : contour["data"][0]) {
                points.push_back({point[0], point[1]});
            }

            // 重新计算长径和短径
            auto [calculated_long, calculated_short] = calculate_diameters(points);

            // 输出结果
            std::cout << "Contour sliceId " << contour["sliceId"] << ":\n";
            std::cout << "  Calculated long diameter: " << calculated_long << "\n";
            std::cout << "  Calculated short diameter: " << calculated_short << "\n";

            // 检查是否与原始值匹配
            double original_long = nodule["longDiameter"];
            double original_short = nodule["shortDiameter"];
            std::cout << "  Original long diameter: " << original_long << "\n";
            std::cout << "  Original short diameter: " << original_short << "\n";
            std::cout << "  Long diameter match?: " << (std::abs(original_long - calculated_long) < 0.1) << "\n";
            std::cout << "  Short diameter match?: " << (std::abs(original_short - calculated_short) < 0.1) << "\n\n";
        }
    }

    return 0;
}



