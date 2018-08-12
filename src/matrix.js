/*
	project onto X/Y/Z

	Eigenwert
	Eigenvektor
*/

class Matrix {
	constructor(values = [], rows, cols) {
		this.rows = rows || undefined;
		this.cols = cols || undefined;
		this.values = [];

		if (values.length > 0) {
			let flat = values.filter(val => Array.isArray(val)).length === 0;
			if (!flat || (typeof rows !== 'undefined' && typeof cols !== 'undefined')) {
				this.fill(values);
			} else {
				throw new Exception('The shape of the matrix could not be determined because the Matrix class was constructed with a flat data array and without the rows and/or cols parameter. Either construct the Matrix class with a two-dimensional array of data in a well-defined shape, or pass the rows and cols parameter on initialisation.');
			}
		}

	}

	fill(values) {
		let flat = values.filter(val => Array.isArray(val)).length === 0;

		if (flat) {
			for (let y = 0; y < this.rows; y++) {
				let row = [];
				for (let x = 0; x < this.cols; x++) {
					row[x] = values[y * this.cols + x] || 0.0;
				}
				this.values[y] = row;
			}
		} else {
			let sameSize = values.filter(val => val.length !== values[0].length).length === 0;
			if (!sameSize) {
				throw new Exception('The shape of the matrix could not be determined because the rows are not of the same length.');
				return;
			}

			this.rows = values.length;
			this.cols = values[0].length;
			this.values = values;
		}

		return this;
	}

	update(matrix) {
		if (typeof matrix === undefined) return this;

		this.values = matrix.values;
		this.rows = matrix.rows;
		this.cols = matrix.cols;
		return this;
	}



	/******************************************************************************
		Static methods invoked on the class e.g. Matrix.add(matrix1, matrix2)
		and public methods invoked on the instance, updating it
	 ******************************************************************************/

	static isMatrix(...matrices) { return matrices.filter(matrix => !(matrix instanceof Matrix)).length === 0; }

	static clone(matrix) { return JSON.parse(JSON.stringify(matrix)); }
	clone() { return Matrix.clone(this); }

	static sameShape(matrix1, ...matrices) {
		if (!Matrix.isMatrix(matrix1, ...matrices)) return;
		return matrices.filter(matrix => matrix.rows !== matrix1.rows || matrix.cols !== matrix1.cols).length == 0;
	}

	sameShape(...matrices) { return Matrix.sameShape(...matrices); }
	
	static log(matrix) { console.table(matrix.values); }
	log() { Matrix.log(this); }



	static add(...matrices) {
		if (!Matrix.sameShape(...matrices)) return;

		return new Matrix(Array(matrices[0].rows).fill(null).map((_, row) =>
			Array(matrices[0].cols).fill(null).map((__, col) =>
				matrices.reduce((acc, val) => acc + val.values[row][col], 0)
			)
		));
	}

	add(...matrices) { return this.update(Matrix.add(this, ...matrices)); }
	


	static addExact(...matrices) {
		if (!Matrix.sameShape(...matrices)) return;

		return new Matrix(Array(matrices[0].rows).fill(null).map((_, row) =>
			Array(matrices[0].cols).fill(null).map((__, col) =>
				sumExact(...matrices.map(matrix => matrix.values[row][col]))
			)
		));
	}

	addExact(...matrices) { return this.update(Matrix.addExact(this, ...matrices)); }
	


	static subtract(matrix1, ...matrices) {
		if (!Matrix.sameShape(matrix1, ...matrices)) return;

		return new Matrix(Array(matrix1.rows).fill(null).map((_, row) =>
			Array(matrix1.cols).fill(null).map((__, col) =>
				matrices.reduce((acc, val) => acc - val.values[row][col], matrix1.values[row][col])
			)
		));
	}

	subtract(...matrices) { return this.update(Matrix.subtract(this, ...matrices)); }
	


	static subtractExact(matrix1, ...matrices) {
		if (!Matrix.sameShape(matrix1, ...matrices)) return;

		return new Matrix(Array(matrix1.rows).fill(null).map((_, row) =>
			Array(matrix1.cols).fill(null).map((__, col) =>
				sumExact(...matrices.map((matrix, i) => i == 0 ? matrix.values[row][col] : - matrix.values[row][col]))
			)
		));
	}

	subtractExact(...matrices) { return this.update(Matrix.subtractExact(this, ...matrices)); }
	


	static mult(matrix, scalar) {
		if (!Matrix.isMatrix(matrix) || !isNumeric(scalar)) return;
		return new Matrix(Array(matrix.rows).fill(null).map((_, row) =>
			Array(matrix.cols).fill(null).map((__, col) => matrix.values[row][col] * scalar)
		));
	}

	mult(scalar) { return this.update(Matrix.mult(this, scalar)); }
	


	static div(matrix, scalar) {
		if (!Matrix.isMatrix(matrix) || !isNumeric(scalar) || scalar == 0) return;
		
		return new Matrix(Array(matrix.rows).fill(null).map((_, row) =>
			Array(matrix.cols).fill(null).map((__, col) => matrix.values[row][col] / scalar)
		));
	}

	div(scalar) { return this.update(Matrix.div(this, scalar)); }
	


	static negative(matrix) {
		if (!Matrix.isMatrix(matrix)) return;
		return Matrix.mult(matrix, -1);
	}

	negative() { return this.update(Matrix.negative(this)); }
	


	static matMult(matrix1, matrix2) {
		if (!Matrix.isMatrix(matrix1, matrix2) || matrix1.cols !== matrix2.rows) return;

		return new Matrix(Array(matrix1.rows).fill(null).map((_, row) =>
			Array(matrix2.cols).fill(null).map((__, col) =>
				Array(matrix1.cols).fill(null).map((__, i) =>
					matrix1.values[row][i] * matrix2.values[i][col]
				).reduce((acc, val) => acc + val, 0)
			)
		))
	}

	matMult(matrix) { return this.update(Matrix.matMult(this, matrix)); }
	


	static submatrix(matrix, rows, cols) {
		if (!Matrix.isMatrix(matrix)) return;
		let newMatrix = Matrix.clone(matrix),
			index;

		newMatrix.values.splice(rows, 1);
		index = newMatrix.values.length;

		while (index--)
			newMatrix.values[index].splice(cols, 1);

		newMatrix.rows = newMatrix.cols = newMatrix.values.length;
		return newMatrix;
	}

	submatrix(rows, cols) { return Matrix.submatrix(this, rows, cols); }
	


	static determinant(matrix) {
		if (!Matrix.isMatrix(matrix) || matrix.rows !== matrix.cols) return;
		if (matrix.rows == 0) return 0;
		if (matrix.rows == 1) return matrix.values[0][0];
		if (matrix.rows == 2) return matrix.values[0][0] * matrix.values[1][1] - matrix.values[0][1] * matrix.values[1][0];

		let sum = 0;
		for (let i = 0; i < 1; i++) {
			for (let j = 0; j < matrix.rows; j++) {
				let sign = (i + j) % 2 ? -1 : 1;
				sum += sign * matrix.values[i][j] * Matrix.determinant(Matrix.submatrix(matrix, i, j));
			}
		}

		return sum;
	}

	determinant() { return Matrix.determinant(this); }


}


