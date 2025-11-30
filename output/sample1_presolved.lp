\ Model _pre
\ LP format - for model browsing. Use MPS format to capture full model detail.
\ Signature: 0xb7ec7bbd9a324be5
Minimize
  y
Subject To
 c1: x1 - 2 r1 <= 3
 c2: - y + 2 r1 + y1 <= -6
 c3: x2 - 3 r2 <= 3
 c4: - y + 3 r2 + y2 <= -5
 c5: x3 - 9 r3 <= 0
 c6: - y + 9 r3 + y3 <= -2
 c7: x1 - 2 r1 - x2 - 13 p12 - 13 q12 <= -8
 c8: 2 r1 + y1 - y2 - 27 p12 + 27 q12 <= 21
 c9: - x1 + x2 - 3 r2 + 14 p12 - 14 q12 <= 6
 c10: - y1 + 3 r2 + y2 + 27 p12 + 27 q12 <= 49
 c11: x1 - 2 r1 - x3 - 13 p13 - 13 q13 <= -8
 c12: 2 r1 + y1 - y3 - 27 p13 + 27 q13 <= 21
 c13: - x1 + x3 - 9 r3 + 11 p13 - 11 q13 <= 0
 c14: - y1 + 9 r3 + y3 + 27 p13 + 27 q13 <= 52
 c15: x2 - 3 r2 - x3 - 14 p23 - 14 q23 <= -8
 c16: 3 r2 + y2 - y3 - 27 p23 + 27 q23 <= 22
 c17: - x2 + x3 - 9 r3 + 11 p23 - 11 q23 <= 0
 c18: - y2 + 9 r3 + y3 + 27 p23 + 27 q23 <= 52
Bounds
 y >= 6
 x1 <= 5
 x2 <= 6
 x3 <= 9
Binaries
 r1 r2 r3 p12 q12 p13 q13 p23 q23
End
