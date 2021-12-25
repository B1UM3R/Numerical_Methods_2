from copy import deepcopy


def Gauss_method(matrix, ans):
    for i in range(len(matrix)):
        a = matrix[i][i]
        for k in range(i, len(matrix)):
            matrix[i][k] /= a
        ans[i] /= a
        for k in range(i + 1, len(matrix)):
            c = -matrix[k][i]
            for j in range(i, len(matrix)):
                matrix[k][j] += c * matrix[i][j]
            ans[k] += c * ans[i]

    for i in range(len(matrix) - 1, -1, -1):
        for k in range(i - 1, -1, -1):
            c = -matrix[k][i]
            matrix[k][i] += c * matrix[i][i]
            ans[k] += c * ans[i]

    return ans.copy()


def Gauss_pivot_selection(matrix, ans):
    n = len(matrix)
    visited = [False] * n
    order = []
    for _ in range(n):
        max_items = map(lambda x: (x[0], *max(enumerate(x[1]), key=lambda y: abs(y[1]))), enumerate(matrix))
        i, j = max(filter(lambda y: not visited[y[0]], max_items), key=lambda x: abs(x[2]))[:2]
        a = matrix[i][j]
        visited[i] = True
        order.append((i, j))
        for k in range(n):
            matrix[i][k] /= a
        ans[i] /= a
        for k in range(n):
            if visited[k]:
                continue
            c = -matrix[k][j]
            for l in range(n):
                matrix[k][l] += c * matrix[i][l]
            ans[k] += c * ans[i]

    answer = [None] * n
    for i in range(n - 1, -1, -1):
        for k in range(i - 1, -1, -1):
            c = -matrix[order[k][0]][order[i][1]]
            matrix[order[k][0]][order[i][1]] += c * matrix[order[i][0]][order[i][1]]
            ans[order[k][0]] += c * ans[order[i][0]]
        answer[order[i][1]] = ans[order[i][0]]

    return answer


def main():
    matrix = [[-0.1, 0.81, 0.38],
              [-0.11, -0.39, -0.5],
              [1.69, -0.01, 0.2]]
    ans = [2.66, -2.39, 2.27]
    result_first = Gauss_method(deepcopy(matrix), deepcopy(ans))
    result_second = Gauss_pivot_selection(deepcopy(matrix), deepcopy(ans))
    print("Решение системы методом Гаусса:", result_first)
    print("Решение системы методом Гаусса c выбором главного элемента по всей матрице:", result_second)


if __name__ == "__main__":
    main()
