Minimize
  20 y + max_n0x - min_n0x + max_n0y - min_n0y
Subject To
  c1: x1 - 2 r1 <= 3
  c2: y1 + 2 r1 - y <= -6
  c3: x2 - 3 r2 <= 3
  c4: y2 + 3 r2 - y <= -5
  c5: x3 - 9 r3 <= 0
  c6: y3 + 9 r3 - y <= -2
  c7: x1 - x2 - 2 r1 - 27 p12 - 27 q12 <= -8
  c8: y1 - y2 + 2 r1 - 27 p12 + 27 q12 <= 21
  c9: x1 - x2 + 3 r2 - 27 p12 + 27 q12 >= -19
  c10: y1 - y2 - 3 r2 - 27 p12 - 27 q12 >= -49
  c11: x1 - x3 - 2 r1 - 27 p13 - 27 q13 <= -8
  c12: y1 - y3 + 2 r1 - 27 p13 + 27 q13 <= 21
  c13: x1 - x3 + 9 r3 - 27 p13 + 27 q13 >= -16
  c14: y1 - y3 - 9 r3 - 27 p13 - 27 q13 >= -52
  c15: x2 - x3 - 3 r2 - 27 p23 - 27 q23 <= -8
  c16: y2 - y3 + 3 r2 - 27 p23 + 27 q23 <= 22
  c17: x2 - x3 + 9 r3 - 27 p23 + 27 q23 >= -16
  c18: y2 - y3 - 9 r3 - 27 p23 - 27 q23 >= -52
  c19: max_n0x - x1 >= 4
  c20: min_n0x - x1 <= 4
  c21: max_n0y - y1 >= 3
  c22: min_n0y - y1 <= 3
  c23: max_n0x - x2 >= 4
  c24: min_n0x - x2 <= 4
  c25: max_n0y - y2 >= 2
  c26: min_n0y - y2 <= 2
  c27: max_n0x - x3 >= 5
  c28: min_n0x - x3 <= 5
  c29: max_n0y - y3 >= 1
  c30: min_n0y - y3 <= 1
  c31: max_n0x >= 0
  c32: min_n0x <= 0
  c33: max_n0y - 0.5 y >= 0
  c34: min_n0y - 0.5 y <= 0
  c35: max_n0x >= 5
  c36: min_n0x <= 5
  c37: max_n0y - y >= 0
  c38: min_n0y - y <= 0
  c39: max_n0x >= 11
  c40: min_n0x <= 11
  c41: max_n0y - 0.5 y >= 0
  c42: min_n0y - 0.5 y <= 0
  c43: max_n0x >= 5
  c44: min_n0x <= 5
  c45: max_n0y >= 0
  c46: min_n0y <= 0
Bounds

Integers

Binaries
  r1 r2 r3 p12 p13 p23 q12 q13 q23
End
