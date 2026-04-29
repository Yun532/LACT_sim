#pragma once

// 分段主镜的“母面”类型
// 注意：这里说的是镜片中心和法向所依附的整体主镜外形，
// 不是镜片本身一定是什么精确曲面。
enum class DishType {
    DaviesCotton,  // DC：镜片中心落在球面母面上
    Parabolic      // 抛物面：镜片中心落在抛物面母面上
};
