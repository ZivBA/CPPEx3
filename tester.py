import numpy
import random
import subprocess
from subprocess import PIPE, STDOUT

PROCESS_NAME = "a.exe"
INTS = 1
DOUBLE = 2
COMPLEX = 3

SUM = 1
MULTI = 2
TRANS = 3


def generateMatrix(rows, cols, comp=INTS):
    mat = []
    for i in range(rows):
        mat.append([])
        for k in range(cols):
            if(comp == COMPLEX):
                x = random.randrange(-10,10)
                y = random.randrange(-10,10)
                mat[i].append(complex(x,y))
                continue
            elif(comp == DOUBLE):
                x = random.randrange(-10,10)
                x -= (x/10)
                mat[i].append(x)
                continue
            mat[i].append(random.randrange(-10,10))

    return numpy.asarray(mat)
def parseToComplex(s):
    minus = False
    sign = s.find("+")
    if(sign == -1):
        minus = True
        sign = s.find("-",1)
    real = int(s[:sign])
    im = int(s[sign+2:-1])
    if(minus):
        im *= -1
    return real, im
def fromOutToMatrix(out,x,y):
    mat = []
    k = False
    relevantLines = []
    for idx,line in enumerate(out.split("\n")):
        if(line == "Resulted matrix:"):
            relevantLines = out.split("\n")[idx+1:idx+x+1]
            break
    for line in relevantLines:
        mat.append(line.split("\t")[:-1])
    if("i" in mat[0][0]):
        k = True
    for line in range(len(mat)):
        for elem in range(len(mat[0])):
            if(k):
                r, i = parseToComplex(mat[line][elem])
                mat[line][elem] = complex(r,i)
                continue
            mat[line][elem] = float(mat[line][elem])

    return numpy.asarray(mat)

def generateToSend(type, oper, x1, y1, matrix1, x2="", y2="", matrix2=numpy.array([])):
    s = ""
    s += str(type) + "\n"
    s += str(oper) + "\n"
    s += str(x1) + "\n" + str(y1) + "\n"
    for list in matrix1:
        for cell in list:
            if (type == COMPLEX):
                if(cell.imag < 0):
                    s += str(cell.real) + "-" + str(cell.imag*-1) + "i,"
                    continue
                s += str(cell.real) + "+" + str(cell.imag) + "i,"
                continue
            s += str(cell)+","
        s = s[:-1]
        s += "\n"
    if(x2 == ""):
        return s
    s += str(x2) + "\n" + str(y2) + "\n"
    for list in matrix2:
        for cell in list:
            if (type == COMPLEX):
                if (cell.imag < 0):
                    s += str(cell.real) + "-" + str(cell.imag * -1) + "i,"
                    continue
                s += str(cell.real) + "+" + str(cell.imag) + "i,"
                continue
            s += str(cell) + "i,"
        s = s[:-1]
        s += "\n"
    return s



def generateMulticipationTest(comp = INTS):
    for i in range(1):
        x = random.randint(1, 20)
        y = random.randint(1, 20)
        z = random.randint(1, 20)
        matrix1 = generateMatrix(x, y, comp)
        matrix2 = generateMatrix(y, z, comp)
        toSend = generateToSend(comp, MULTI, x, y, matrix1, y, z, matrix2)
        p = subprocess.Popen([PROCESS_NAME], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        output, error = p.communicate(input=toSend.encode('utf-8'), timeout=15)
        matrix3 = fromOutToMatrix(output.decode('utf-8'), x, z)
        matrix4 = numpy.dot(matrix1, matrix2)
        if (not numpy.allclose(matrix3, matrix4)):
            print("error by trying to MUL")
            print(matrix1)
            print(matrix2)
            print("user output")
            print(matrix3)
            print("numpy output")
            print(matrix4)

def generateSummingTest(comp = INTS):
    for i in range(1):
        x = random.randint(1, 3)
        y = random.randint(1, 3)
        matrix1 = generateMatrix(x, y, comp)
        matrix2 = generateMatrix(x, y, comp)
        toSend = generateToSend(comp, SUM, x, y, matrix1, x, y, matrix2)
        p = subprocess.Popen([PROCESS_NAME], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        output, error = p.communicate(input=toSend.encode('utf-8'), timeout=15)
        matrix3 = fromOutToMatrix(output.decode('utf-8'), x, y)
        if (not numpy.allclose(matrix3, matrix1+matrix2)):
            print(toSend)
            print("OUTPUT")
            print(output.decode())
            print("error by trying to SUM")
            print(matrix1)
            print(matrix2)
            print("user output")
            print(matrix3)
            print("real output")
            print(matrix1 + matrix2)

def generateTransCheck(comp = INTS):
    for i in range(1):
        x = random.randint(1, 35)
        y = random.randint(1, 35)
        matrix1 = generateMatrix(x, y, comp)
        toSend = generateToSend(comp, TRANS, x, y, matrix1)
        p = subprocess.Popen([PROCESS_NAME], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        output, error = p.communicate(input=toSend.encode('utf-8'), timeout=15)
        matrix3 = fromOutToMatrix(output.decode('utf-8'), y,x)
        if (not numpy.array_equal(matrix3, matrix1.transpose().conj())):
            print(output.decode())
            print("error by trying to TRANS")
            print(matrix1)
            print("user output")
            print(matrix3)
            print("numpy output")
            print(matrix1.transpose().conj())
        return matrix3

for i in range(100):
    print("Loop number "+str(i))
    generateTransCheck(INTS)
    generateTransCheck(DOUBLE)
    generateTransCheck(COMPLEX)
    generateSummingTest(INTS)
    generateSummingTest(DOUBLE)
    generateSummingTest(COMPLEX)
    generateMulticipationTest(INTS)
    generateMulticipationTest(DOUBLE)
    generateMulticipationTest(COMPLEX)
