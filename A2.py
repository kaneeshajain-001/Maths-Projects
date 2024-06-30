import numpy as np
import math

def read_input(filename):
    with open(filename, 'r') as file:
        content = file.read().splitlines()
    

    # Find the indices for sections
    obj_index = content.index("[objective]") + 1
    A_index = content.index("[A]") + 1
    b_index = content.index("[b]") +1
    constraint_types_index = content.index("[constraint_types]") + 1
    c_index = content.index("[c]") + 1
    #print(content[b_index:constraint_types_index-2])
    # Extract objective, A, b, constraint_types, and c
    objective = content[obj_index].strip()
    
    # Filter out empty strings before splitting
    A = [list(map(float, filter(None, line.split(',')))) for line in content[A_index:b_index-2]]

    b = list(map(float, content[b_index:constraint_types_index-2]))
    constraint_types = content[constraint_types_index:c_index-2]

    
    # Filter out empty strings before splitting
    c = list(map(float, filter(None, content[c_index].split(','))))
    return np.array(A), np.array(b), np.array(c), constraint_types, objective


def is_integer_solution(solution, tolerance=1e-6):
    for value in solution:
        if abs(value - round(value)) > tolerance:
            return False
    return True

def dual_simplex(gomory_cut,tableau,cut_number,zeroth_row,zeroth_col,A):
    precisions=1e-6
    rows,cols=tableau.shape
    solution_status="optimal"
    new_tableau = np.zeros((rows+1, cols+ 1 ))
    for i in range(0,rows):
        for j in range(0,cols):
            new_tableau[i][j]=tableau[i][j]
    for k in range(0,len(gomory_cut)):
        new_tableau[rows][k]=gomory_cut[k]
    
    #print(new_tableau)
    while np.min(new_tableau[1:, 0]) < (-1)*precisions:
        i=1
        while(i<len(new_tableau) and new_tableau[i][0]>=(-1)*precisions):
            i+=1
        pivot_row=i
        pivot_row_index=zeroth_col[i-1]
        pivot_col=0
        pivot_col_index=""
        mini = float('inf')
        
        for j in range(1, len(new_tableau[pivot_row])):
            if new_tableau[pivot_row][j] <(-1)*precisions:
                ratio = (new_tableau[0][j]) / (abs(new_tableau[pivot_row][j]))
                if ratio < mini:
                    mini = ratio
                    pivot_col = j
                    
                    pivot_col_index = zeroth_row[j - 1]
        pivot_element = new_tableau[pivot_row][pivot_col]
     
        for i in range(len(new_tableau[0])):
            new_tableau[pivot_row][i] /= pivot_element
        for i in range(len(new_tableau)):
            if i != pivot_row:
                multiplier = new_tableau[i][pivot_col]
                for j in range(len(new_tableau[0])):
                    new_tableau[i][j] -= multiplier * new_tableau[pivot_row][j]
        
        zeroth_col[pivot_row - 1] = pivot_col_index
    #print(new_tableau)
    sol=np.zeros(len(A[0]))
    #print(zeroth_col)
    for i in range(len(zeroth_col)):
        if zeroth_col[i].startswith('x'):
            sol[int(zeroth_col[i][-1])] =new_tableau[i+1][0]
    #print(sol)
    is_integer=is_integer_solution(sol)
    #print(new_tableau)
    answer = (
        new_tableau, is_integer,zeroth_row,zeroth_col,solution_status
    )

    return answer

def simplex_algo(A, b, c, constraint_types, objective):
    solution_status = ""
    if objective == "maximize":
        for i in range(len(c)):
            c[i] *= -1
    for i in range(len(b)):
        if b[i] < 0:
            if constraint_types[i] == "<=":
                constraint_types[i] = ">="
            elif constraint_types[i] == "<":
                constraint_types[i] = ">"
            elif constraint_types[i] == ">":
                constraint_types[i] = "<"
            elif constraint_types[i] == ">=":
                constraint_types[i] = "<="
            A[i] = -1 * A[i]
            b[i] = -1 * b[i]
    non_artificialvar = 0
    initial_bfs = np.array([0 * (len(A))])
    zeroth_row = np.array([""] * len(A[0]))
    zeroth_row = ["x" + str(i) for i in range(len(zeroth_row))]
    non_artificialvar = non_artificialvar + len(zeroth_row)
    zeroth_col = np.empty(0)
    initial_cost = 0

    cb = np.zeros(len(constraint_types))
    for i in range(len(A)):
        if constraint_types[i] == "<":
            zeroth_row = np.concatenate((zeroth_row, [f"s{str(i)}"]))
            zeroth_col = np.concatenate((zeroth_col, [f"s{str(i)}"]))
            non_artificialvar = non_artificialvar + 1
            cb[i] = 0

        elif constraint_types[i] == "<=":
            zeroth_row = np.concatenate((zeroth_row, [f"s{str(i)}"]))
            zeroth_col = np.concatenate((zeroth_col, [f"s{str(i)}"]))
            non_artificialvar = non_artificialvar + 1
            cb[i] = 0
        elif constraint_types[i] == ">":
            zeroth_row = np.concatenate((zeroth_row, [f"s{str(i)}"]))
            zeroth_row = np.concatenate((zeroth_row, [f"A{str(i)}"]))
            zeroth_col = np.concatenate((zeroth_col, [f"A{str(i)}"]))
            non_artificialvar = non_artificialvar + 1
            cb[i] = 1
            initial_cost = initial_cost + b[i]
        elif constraint_types[i] == ">=":
            zeroth_row = np.concatenate((zeroth_row, [f"s{str(i)}"]))
            zeroth_row = np.concatenate((zeroth_row, [f"A{str(i)}"]))
            zeroth_col = np.concatenate((zeroth_col, [f"A{str(i)}"]))
            non_artificialvar = non_artificialvar + 1
            cb[i] = 1
            initial_cost = initial_cost + b[i]
        elif constraint_types[i] == "=":
            zeroth_row = np.concatenate((zeroth_row, [f"A{str(i)}"]))
            zeroth_col = np.concatenate((zeroth_col, [f"A{str(i)}"]))
            cb[i] = 1
            initial_cost = initial_cost + b[i]

    phase1 = np.zeros((len(zeroth_col) + 1, len(zeroth_row) + 1))
    m, n = phase1.shape
    phase1[0][0] = -1 * initial_cost
    for i in range(1, m):
        phase1[i][0] = b[i - 1]
    for i in range(1,len(A[0])+1):
        cost = 0
        for j in range(1,len(phase1)):
            phase1[j][i] = A[j-1][i-1]
            cost += cb[j-1]*A[j-1][i-1]
        phase1[0][i] = -1*cost
    col_no = len(A[0]) + 1
    for i in range(len(A)):
        if constraint_types[i] == "<":
            phase1[i + 1][col_no] = 1
            phase1[0][col_no]=0
            col_no += 1
        elif constraint_types[i] == "<=":
            phase1[i + 1][col_no] = 1
            phase1[0][col_no]=0
            col_no += 1
        elif constraint_types[i] == ">":
            phase1[i + 1][col_no] = -1
            phase1[i + 1][col_no + 1] = 1
            phase1[0][col_no]=1
            phase1[0][col_no+1]=0
            col_no += 2
        elif constraint_types[i] == ">=":
            phase1[0][col_no] = cb[i] * 1
            phase1[i + 1][col_no] = -1
            phase1[i + 1][col_no + 1] = 1
            phase1[0][col_no]=1
            phase1[0][col_no+1]=0
            col_no += 2
        elif constraint_types[i] == "=":
            phase1[i + 1][col_no] = 1
            phase1[0][col_no]=0
            col_no += 1
    np.set_printoptions(precision=2, suppress=True)
    initial_tableau = phase1

    tolerance = 1e-9
    while any(phase1[0, 1:] < 0) and abs(phase1[0][0]) > tolerance:
        pivot_col = 0
        pivot_col_index = ""
        for i in range(1, len(phase1[0])):
            if phase1[0][i] < 0:
                pivot_col = i
                pivot_col_index = zeroth_row[i - 1]
                break
        pivot_row = 0
        pivot_row_index = ""
        mini = float('inf')
        for i in range(1, len(phase1)):
            if phase1[i][pivot_col] > 0:
                ratio = phase1[i][0] / phase1[i][pivot_col]
                if ratio < mini:
                    mini = ratio
                    pivot_row = i
                    pivot_row_index = zeroth_col[i - 1]
        pivot_element = phase1[pivot_row][pivot_col]
        for i in range(len(phase1[0])):
            phase1[pivot_row][i] /= pivot_element
        for i in range(len(phase1)):
            if i != pivot_row:
                multiplier = phase1[i][pivot_col]
                for j in range(len(phase1[0])):
                    phase1[i][j] -= multiplier * phase1[pivot_row][j]

        zeroth_col[pivot_row - 1] = pivot_col_index
        # idhar wala tableau final hai

    tolerance = 1e-9

    if abs(phase1[0][0])> tolerance:
        solution_status = "infeasible"
    else:
        while any(col.startswith('A') for col in zeroth_col):
            i = 0
            while(i < len(zeroth_col)):
                if zeroth_col[i].startswith('A'):
                    pivot_row = i + 1
                    pivot_row_index = zeroth_col[i]
                    flag = True
                    for j in range(1, non_artificialvar + 1):
                        if phase1[pivot_row][j] != 0:
                            flag = False
                            pivot_col = j
                            pivot_col_index = zeroth_row[j - 1]
                            break
                    if not flag:
                        pivot_element = phase1[pivot_row][pivot_col]
                        for i in range(len(phase1[0])):
                            phase1[pivot_row][i] /= pivot_element

                        for i in range(len(phase1)):
                            if i != pivot_row:
                                multiplier = phase1[i][pivot_col]
                                for j in range(len(phase1[0])):
                                    phase1[i][j] -= multiplier * phase1[pivot_row][j]
                        zeroth_col[pivot_row - 1] = pivot_col_index
                    else:
                        phase1 = np.delete(phase1, pivot_row, axis=0)
                        zeroth_col = np.delete(zeroth_col, pivot_row - 1, axis=0)

                else:
                    i = i + 1

        while any(col.startswith('A') for col in zeroth_row):
            for i in range(len(zeroth_row)):
                if zeroth_row[i].startswith('A'):
                    phase1 = np.delete(phase1, i + 1, axis=1)
                    zeroth_row = np.delete(zeroth_row, i)
                    break
    m = len(phase1)
    n = len(phase1[0])
    initial_cost2 = 0
    for i in range(len(zeroth_col)):
        if zeroth_col[i].startswith('x'):
            initial_cost2 = c[int(zeroth_col[i][-1])] * phase1[i + 1][0] + initial_cost2
    phase1[0][0] = -1 * initial_cost2

    basic_cost_vector = [0] * len(zeroth_col)
    for i in range(len(zeroth_col)):
        if zeroth_col[i].startswith('x'):
            basic_cost_vector[i] = c[int(zeroth_col[i][-1])]
        else:
            basic_cost_vector[i] = 0
    for i in range(0, len(zeroth_row)):
        if zeroth_row[i] in zeroth_col:
            phase1[0][i + 1] = 0
        else:
            if zeroth_row[i].startswith('x'):
                cost = c[int(zeroth_row[i][-1])]
                astha = 0
                for j in range(1, len(phase1)):
                    astha += basic_cost_vector[j - 1] * phase1[j][i + 1]
                phase1[0][i + 1] = cost - astha
            else:
                cost = 0
                astha = 0
                for j in range(1, len(phase1)):
                    astha += basic_cost_vector[j - 1] * phase1[j][i + 1]
                phase1[0][i + 1] = cost - astha

    while any(phase1[0, 1:] < 0):
        pivot_col = 0
        pivot_col_index = ""
        for i in range(1, len(phase1[0])):
            if phase1[0][i] < 0:
                pivot_col = i
                pivot_col_index = zeroth_row[i - 1]
                break

        pivot_row = 0
        pivot_row_index = ""
        mini = float('inf')
        found_pivot = False
        for i in range(1, len(phase1)):
            if phase1[i][pivot_col] > 0:
                found_pivot = True
                ratio = phase1[i][0] / phase1[i][pivot_col]
                if ratio < mini:
                    mini = ratio
                    pivot_row = i
                    pivot_row_index = zeroth_col[i - 1]
        if found_pivot == False:
            solution_status = "Unbounded"
            break
        pivot_element = phase1[pivot_row][pivot_col]
        for i in range(len(phase1[0])):
            phase1[pivot_row][i] /= pivot_element
        for i in range(len(phase1)):
            if i != pivot_row:
                multiplier = phase1[i][pivot_col]
                for j in range(len(phase1[0])):
                    phase1[i][j] -= multiplier * phase1[pivot_row][j]

        zeroth_col[pivot_row - 1] = pivot_col_index
    optimal_solution = [0] * len(A[0])
    if solution_status != "Unbounded":
        for i in range(1, len(phase1)):
            if zeroth_col[i - 1].startswith('x'):
                solution_status = "feasible"
                optimal_solution[int(zeroth_col[i - 1][-1])] = phase1[i][0]
    final_tableau = phase1
    if objective == "maximize":
        q = phase1[0][0]
    else:
        q = -1 * phase1[0][0]
    answer = (
        final_tableau, solution_status, optimal_solution,q,zeroth_row,zeroth_col
    )

    return answer


def gomory_cut_algo():
    filename=r'input_ilp.txt' 
    A,b,c,constraint,objective=read_input(filename)
    # A=np.array([[1,-1],[2,4]])
    # b=np.array([2,15])
    # c=np.array([1,-3])    
    # constraint=["<=","<="]
    # objective="minimize"
    cost=np.copy(c)
    ft,ss,os,ov,zeroth_row,zeroth_col=simplex_algo(A,b,c,constraint,objective) # assign variables, ft=final tableau, ss=soln status, os=final optimal soln, ov=opt val
    length=A.shape[0]
    if (ss=="infeasible"):
        solution_status="infeasible"
        return solution_status
    elif (ss=="unbounded"):
        solution_status="unbounded"
        return solution_status
    else:
        initial_solution=os
        isint=True
        for i in range(1,len(A)+1):
            if(ft[i,0]%1<=1e-6 or 1-ft[i,0]%1<=1e-6):
                continue
            else:
                isint=False
                break

        if (isint):
            optimal_value=ov
        else:
            number_of_cuts=0
            source_row=0
            while(isint==False):
                maxi=0;             
                for i in range(1,len(A)+1):
                    if ft[i,0]%1>1e-6:
                        if(ft[i,0]>=maxi):
                    
                            maxi=ft[i,0]

                            source_row=ft[i,:]
                            break
                source_row_size = len(source_row)
                gomory_cut = np.zeros(source_row_size + 1)
                next_index = 'g' + str(number_of_cuts)
                zeroth_row = np.append(zeroth_row, next_index)
                zeroth_col=np.append(zeroth_col,next_index)

                gomory_cut[-1]=1
                gomory_cut[0]=-1*(source_row[0]%1)
                for j in range(1, len(gomory_cut)-1):
                    if zeroth_row[j-1] not in zeroth_col:
                        if(ft[i][j]%1 > (1e-6) and 1-ft[i][j]%1 > (1e-6)):
                            gomory_cut[j] = (-1) * (source_row[j] % 1)
                number_of_cuts=number_of_cuts+1
                ft,isint,zeroth_row, zeroth_col,solution_status=dual_simplex(gomory_cut,ft,number_of_cuts,zeroth_row,zeroth_col,A)
 
            if(isint):
                
                #print(zeroth_col)
                optimal_solution = [0] * len(A[0])
                for i in range(1, len(ft)):
                    if zeroth_col[i - 1].startswith('x'):
                        optimal_solution[int(zeroth_col[i - 1][-1])] = np.round(ft[i][0])
                # optimal value
                optimal_value=0
                
                for i in range(len(cost)):
                    optimal_value+=cost[i]*optimal_solution[i]
                print("initial_solution: ",initial_solution)
                print("final_solution: ",optimal_solution)
                print("solution_status: ",solution_status)
                print("number_og_cuts: ",number_of_cuts)
                print("optimal_value: ",np.round(optimal_value))

        #If the optimal solution is integers, then problem is solved. Otherwise, 
        #add Gomory's constraint (cut) is added to optimal solution. Now new problem 
        #is solved using dual simplex method The method terminates as soon as optimal solution become integers.

        
    # check when the ss is unbdd, infeasible
gomory_cut_algo()
