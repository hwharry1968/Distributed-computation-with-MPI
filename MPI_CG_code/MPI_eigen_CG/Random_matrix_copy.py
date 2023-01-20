import random
import sys

## Generates a random dense matrix
## Use: python randomMatrix.py M N

if (len(sys.argv) > 1):

  M = int(sys.argv[1])
  N = M

  ans = [ 0 for i in range(M) ]
  b = [ 0 for i in range(M) ]
  table2 = [ [ 0 for i in range(N) ] for j in range(M) ]

  for i in range(M):
    ans[i] = random.random()
    b[i] = 0
  
  for i in range(M):
    for j in range(N):
      if i <= j:
        table2[i][j] = random.random()
      if i < j:
        if random.random() < 0.9:
          table2[i][j] = 0
        table2[j][i] = table2[i][j]


  nnz = 0 
  val = []
  col_ind = []
  row_ptr = [0]
  for i in range(M):
    cnt_row_nnz = 0
    for j in range(N):
      if(table2[i][j] != 0):
        b[i] += table2[i][j] * ans[j]
        cnt_row_nnz += 1
        val.append(table2[i][j])
        col_ind.append(j)
    row_ptr.append(cnt_row_nnz)
    nnz += cnt_row_nnz
  
  for i in range(M):
    row_ptr[i + 1] += row_ptr[i]

  f = open('csrmatrix20','w')
  f.write(str(M))
  f.write('\t')
  f.write(str(N))
  f.write('\t')
  f.write(str(nnz))
  f.write('\n')

  for i in row_ptr:
    f.write(str(i))
    f.write('\t')
  f.write('\n')

  for i in col_ind:
    f.write(str(i))
    f.write('\t')
  f.write('\n')

  for i in val:
    f.write(str(i))
    f.write('\t')   
  f.write('\n')

  for i in range(M):
    f.write(str(b[i]))
    f.write('\t')
  f.write('\n')
  ##f.close()
  
  ##f = open('ans','w')
  for i in range(M):
    f.write(str(ans[i]))
    f.write('\t')
  f.write('\n')
  f.close()


else:
  print("please input a number")

