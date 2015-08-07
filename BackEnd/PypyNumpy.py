class vector:
    ''' Vector of unknown length. Inputs is a list. '''
    def __init__(self,r):
        self.r = r
    def __len__(self):
        ''' Returns length of vector. '''
        '''Ex: len(Vector([1,2])) '''
        return (len(self.r))
    def __getitem__(self,key):
        ''' Returns element in vector. '''
        '''Ex: Vector([1,2])[0] '''
        return self.r[key]
    def __setitem__(self,key,value):
        ''' Sets element in vector to a value.'''
        '''Ex: Vector([1,2])[0]=3'''
        if key < len(self.r):
            self.r[key]=value
        if key == len(self.r):
            self.r.append(value)
        else:
            raise IndexError('Invalid index key')
    def __repr__(self):
        ''' Returns a string version of vector.'''
        '''Ex: x=Vector([1,2]) ... x---> Vector [1,2]'''
        return 'vector %s' % self.r
    def __add__(self,other):
        ''' Returns the addition, allows other to be vector or a int/float'''
        ''' Ex: (x is vectors, y is vector/int/float) z = x+y '''
        if isinstance(other,vector) or isinstance(other,list):
            return vector([self.r[i] + other[i] for i in range(len(self.r))])
        if isinstance(other,float) or isinstance(other,int):
            return vector([self.r[i] + other for i in range(len(self.r))])
        else:
            raise TypeError('other must be Vector,float, or int')
    def __radd__(self,other):
        ''' Returns the addition, allows other to be vector or a int/float'''
        ''' Ex: (x is vectors, y is vector/int/float) z = x+y '''
        if isinstance(other,vector) or isinstance(other,list):
            return vector([self.r[i] + other[i] for i in range(len(self.r))])
        if isinstance(other,float) or isinstance(other,int):
            return vector([self.r[i] + other for i in range(len(self.r))])
        else:
            raise TypeError('other must be Vector,float, or int')
    def __iadd__(self,other):
        if isinstance(other,vector) or isinstance(other,list):
            return vector([self.r[i] + other[i] for i in range(len(self.r))])
        if isinstance(other,float) or isinstance(other,int):
            return vector([self.r[i] + other for i in range(len(self.r))])
    def __sub__(self,other):
        ''' Returns the addition, allows other to be vector or a int/float'''
        ''' Ex: (x is vectors, y is vector/int/float) z = x-y '''
        if isinstance(other,vector) or isinstance(other,list):
            return vector([self.r[i] - other[i] for i in range(len(self.r))])
        if isinstance(other,float) or isinstance(other,int):
            return vector([self.r[i] - other for i in range(len(self.r))])
        else:
            raise TypeError('other must be Vector,list,float, or int')
    def __rsub__(self,other):
        ''' Returns the addition, allows other to be vector or a int/float'''
        ''' Ex: (x is vectors, y is vector/int/float) z = x-y '''
        if isinstance(other,vector) or isinstance(other,list):
            return vector([other[i] - self.r[i] for i in range(len(self.r))])
        if isinstance(other,float) or isinstance(other,int):
            return vector([other - self.r[i] for i in range(len(self.r))])
        else:
            raise TypeError('other must be Vector,list,float, or int')
    def __isub__(self,other):
        if isinstance(other,vector) or isinstance(other,list):
            return vector([self.r[i] - other[i] for i in range(len(self.r))])
        if isinstance(other,float) or isinstance(other,int):
            return vector([self.r[i] - other for i in range(len(self.r))])
        else:
            raise TypeError('other must be Vector,list,float, or int')
    def __mul__(self,other):
        ''' Returns the element by element mulutiplication'''
        ''' x*y = z '''
        if isinstance(other,vector) or isinstance(other,list):
            return vector([self.r[i]*other[i] for i in range(len(self.r))])
        if isinstance(other,float) or isinstance(other,int):
            return vector([self.r[i]*other for i in range(len(self.r))])
        else:
            raise TypeError('other must be Vector or list')
    def __rmul__(self,other):
        ''' Returns the element by element mulutiplication'''
        ''' x*y = z '''
        if isinstance(other,vector) or isinstance(other,list):
            return vector([self.r[i]*other[i] for i in range(len(self.r))])
        if isinstance(other,float) or isinstance(other,int):
            return vector([self.r[i]*other for i in range(len(self.r))])
        else:
            raise TypeError('other must be Vector or list')
    def __imul__(self,other):
        if isinstance(other,vector) or isinstance(other,list):
            return vector([self.r[i]*other[i] for i in range(len(self.r))])
        if isinstance(other,float) or isinstance(other,int):
            return vector([self.r[i]*other for i in range(len(self.r))])
        else:
            raise TypeError('other must be Vector or list')
    def __truediv__(self,other):
        ''' Returns the element by element division'''
        ''' x.divide(y) = z '''
        if isinstance(other,vector) or isinstance(other,list):
            return vector([self.r[i]/other[i] for i in range(len(self.r))])
        if isinstance(other,float) or isinstance(other,int):
            return vector([self.r[i]/other for i in range(len(self.r))])
        else:
            raise TypeError('other must be Vector or list')
    def __rtruediv__(self,other):
        ''' Returns the element by element division'''
        ''' x.divide(y) = z '''
        if isinstance(other,vector) or isinstance(other,list):
            return vector([other[i]/self.r[i] for i in range(len(self.r))])
        if isinstance(other,float) or isinstance(other,int):
            return vector([other/self.r[i] for i in range(len(self.r))])
        else:
            raise TypeError('other must be Vector or list')
    def __itruediv__(self,other):
        if isinstance(other,vector) or isinstance(other,list):
            return vector([self.r[i]/other[i] for i in range(len(self.r))])
        if isinstance(other,float) or isinstance(other,int):
            return vector([self.r[i]/other for i in range(len(self.r))])
        else:
            raise TypeError('other must be Vector or list')
    def __pow__(self,power):
        ''' Returns the vector raised to a power '''
        ''' Ex: x.pow(2) '''
        assert type(power) == int or type(power) == float, 'Must be raised to a int or float'
        return vector([self.r[i]**power for i in range(len(self.r))])
    def average(self):
        ''' Returns average of vector elements'''
        ''' Ex: x.averager() '''
        return sum(self.r)/len(self.r)
    def norm(self):
        ''' Returns norm of vector'''
        ''' Ex: x.norm() '''
        return (sum([self.r[i]*self.r[i] for i in range(len(self.r))]))**(1/2)
    def dot(self,other):
        ''' Returns dot product of self and other '''
        assert len(other) == len(self.r) , 'Vector lengths must be equal'
        return sum([a*b for a,b in zip(self.r,other)])
    def exp(self):
        e = 2.718281828459045
        return vector([e**a for a in self.r])

class matrix:
    def __init__(self,x):
        self.x = x
    def __str__(self):
        s = ""
        for row in self.x:
            s += "%s\n" % row
        return s
    def __getitem__(self,key):
        ''' returns a value in the matrix '''
        ''' Ex: m = matrix([[1,5],[1,1]]) ... m[0,1] = 5 '''
        row , col = key
        return self.x[row][col]
    def __setitem__(self,key,value):
        ''' ammends a value in the matrix '''
        row , col = key
        self.x[row][col] = value
    def __repr__(self):
        return 'matrix %s' % self.x
    def dimension(self):
        ''' returns a list with the number of rows and columns of the matrix '''
        return [len(self.x),len(self.x[0])]
    def transpose(self):
        ''' returns transpose of matrix '''
        return matrix([list(i) for i in zip(*self.x)])
    def __add__(self,other):
        ''' adds two matrix's element by element '''
        assert type(other) == matrix , 'Must add matrix with matrix not matrix'
        assert other.dimension() == [len(self.x),len(self.x[0])], 'Matrix dimensions do not agree'
        return matrix([[a+b for a,b in zip(self.x[row],other[row,:])] for row in range(len(self.x[0]))])
    def __radd__(self,other):
        ''' adds a matrix from the left '''
        assert type(other) == matrix , 'Must add matrix with matrix not matrix'
        assert other.dimension() == [len(self.x),len(self.x[0])], 'Matrix dimensions do not agree'
        return matrix([[a+b for a,b in zip(other[row,:],self.x[row])] for row in range(len(self.x[0]))])
    def __iadd__(self,other):
        ''' adds a matrix to this matrix'''
        assert type(other) == matrix , 'Must add matrix with matrix not matrix'
        assert other.dimension() == [len(self.x),len(self.x[0])], 'Matrix dimensions do not agree'
        return matrix([[a+b for a,b in zip(self.x[row],other[row,:])] for row in range(len(self.x[0]))])
    def __sub__(self,other):
        ''' subtracts two matrix's element by element '''
        assert type(other) == matrix , 'Must subtract matrix with matrix not matrix'
        assert other.dimension() == [len(self.x),len(self.x[0])], 'Matrix dimensions do not agree'
        return matrix([[a-b for a,b in zip(self.x[row],other[row,:])] for row in range(len(self.x[0]))])
    def __rsub__(self,other):
        ''' adds a matrix from the left '''
        assert type(other) == matrix , 'Must subtract matrix with matrix not matrix'
        assert other.dimension() == [len(self.x),len(self.x[0])], 'Matrix dimensions do not agree'
        return matrix([[a-b for a,b in zip(other[row,:],self.x[row])] for row in range(len(self.x[0]))])
    def __isub__(self,other):
        ''' adds a matrix to this matrix'''
        assert type(other) == matrix , 'Must subtract matrix with matrix not matrix'
        assert other.dimension() == [len(self.x),len(self.x[0])], 'Matrix dimensions do not agree'
        return matrix([[a-b for a,b in zip(self.x[row],other[row,:])] for row in range(len(self.x[0]))])
    def __mul__(self,other):
        ''' multiplies two matrixes element by element '''
        if type(other) == int or type(other) == float:
            return matrix([[a*other for a in self.x[X_row] ]for X_row in range(len(self.x))])
        else:
            assert len(self.x[0]) ==  other.dimension()[0] , 'Matrix dimensions do not agree'
            return matrix([[sum(a*b for a,b in zip(X_row,Y_col)) for Y_col in zip(*other[:,:])] for X_row in self.x])
    def __rmul__(self,other):
        ''' multiplies two matrixes element by element from the left '''
        if type(other) == int or type(other) == float:
            return matrix([[a*other for a in self.x[X_row] ]for X_row in range(len(self.x))])
        else:
            assert other.dimension()[0] == len(self.x[:,0]), 'Matrix dimensions do not agree'
            return matrix([[sum(a*b for a,b in zip(X_row,Y_col)) for Y_col in zip(*self.x[:,:])] for X_row in other[:,:]])
    def __imul__(self,other):
        ''' multiplies a matrix to this matrix'''
        if type(other) == int or type(other) == float:
            return matrix([[a*other for a in self.x[X_row] ]for X_row in range(len(self.x))])
        else:
            assert other.dimension() == [len(self.x),len(self.x[0])], 'Matrix dimensions do not agree'
            return matrix([[sum(a*b for a,b in zip(X_row,Y_col)) for Y_col in zip(*other[:,:])] for X_row in self.x])
    def __truediv__(self,other):
        ''' Divides to matrixs '''
        if type(other) == int or type(other) == float:
            return matrix([[a/other for a in self.x[X_row] ]for X_row in range(len(self.x))])
        else:
            assert len(self.x[0]) ==  other.dimension()[0] , 'Matrix dimensions do not agree'
            return matrix([[a/b for a,b in zip(self.x[row],other[row,:])] for row in range(len(self.x[0]))])
    def __rtruediv__(self,other):
        ''' Divides to matrixs '''
        if type(other) == int or type(other) == float:
            return matrix([[other/a for a in self.x[X_row] ]for X_row in range(len(self.x))])
        else:
            assert other.dimension()[0] == len(self.x[:,0]), 'Matrix dimensions do not agree'
            return matrix([[(a/b for a,b in zip(X_row,Y_col)) for Y_col in zip(*self.x[:,:])] for X_row in other[:,:]])
    def __itruediv__(self,other):
        ''' Divides to matrixs '''
        if type(other) == int or type(other) == float:
            return matrix([[a/other for a in self.x[X_row] ]for X_row in range(len(self.x))])
        else:
            assert other.dimension() == [len(self.x),len(self.x[0])], 'Matrix dimensions do not agree'
            return matrix([[a/b for a,b in zip(self.x[row],other[row,:])] for row in range(len(self.x[0]))])
    def exp(self):
        e = 2.718281828459045
        return matrix([[e**a for a in self.x[X_row] ]for X_row in range(len(self.x))])


