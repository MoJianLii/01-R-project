# 第1步：设置编译环境变量
Sys.setenv("PKG_CXXFLAGS" = "-std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0")
Sys.setenv("PKG_CPPFLAGS" = "-I/usr/include")

# 第2步：重新测试
library(Rcpp)
evalCpp("2 + 2")  # 应该返回4

cppFunction('int add(int x, int y) { return x + y; }')
add(2, 3)  # 应该返回5
# 立即运行Rcpp版本（最快！）
source("sample_500_rcpp_FIXED_PRECISION_svr4.R")
# 预期：3-6分钟/sample，总共1.4天
