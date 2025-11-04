# LGA_hotFIX for Local-Global Alignment (LGA) program
发现用于计算GDT_TS的代码存在栈溢出的问题，在特定设备上会崩溃，随即进行修复。

## 用法
make -f Makefile.linux_x86_64_msse2 clean && make -f Makefile.linux_x86_64_msse2

直接编译就可以，用编译出来的lga执行计算。
