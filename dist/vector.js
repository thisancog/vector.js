var isNumeric = function(n) {
	return (typeof n !== 'undefined') ? !Array.isArray(n) && !isNaN(parseFloat(n)) && isFinite(n) : undefined;
}

var sumExact = function(data) {
	if (data.length == 0) return;
	if (data.length == 1) return data[0];

	let sum = data[0],
		compensation = 0.0,
		decimals = 0;

	for (let i = 1; i < data.length; i++) {
		let summand = data[i],
			t = sum + summand;

		if ((summand + "").indexOf('.') > -1) {
			let numDec = (summand + "").split(".")[1].length;
			decimals = numDec > decimals ? numDec : decimals;
		}

		compensation += (Math.abs(sum) >= Math.abs(summand)) ? (sum - t) + summand : (summand - t) + sum;
		sum = t;
	}

	let power = Math.pow(10, decimals + 1);
	return Math.round((sum + compensation) * power) / power;
}
/*
	matrix multiplication
	matrix division
	rotation matrix X
	rotation matrix Y
	rotation matrix Z
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



	/******************************************
		Static methods invoked on the class
		e.g. Matrix.add(matrix1, matrix2)
	 ******************************************/

	static submatrix(matrix, rows, cols) {
		let newMatrix = Matrix.clone(matrix),
			index;

		newMatrix.values.splice(rows, 1);
		index = newMatrix.values.length;

		while (index--)
			newMatrix.values[index].splice(cols, 1);

		newMatrix.rows = newMatrix.cols = newMatrix.values.length;
		return newMatrix;
	}

	static determinant(matrix) {
		if (matrix.rows !== matrix.cols) return;
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

	static clone(matrix) {
		return JSON.parse(JSON.stringify(matrix));
	}



	/******************************************
		Public methods invoked on the instance
		which update the instance's data
	 ******************************************/

	determinant() {
		return Matrix.determinant(this);
	}

	submatrix(rows, cols) {
		return Matrix.submatrix(this, rows, cols);
	}

	clone() {
		return Matrix.clone(this);
	}

}




/*

	ToDo:

	cartesian to polar
	polar to cartesian
	construct from polar
	direction in cartesian grid
	angle
	angle in plane N

	cross product
	scalar triple product
	vector triple product
	dyadic product

	projection a on b

*/

"use strict";

class Vector {
	constructor(...args) {
		let dim = 0;

		args.forEach(function(arg) {
			if (isNumeric(arg)) {
				let label = Vector.getDimensionLabel(dim);
				this[label] = parseFloat(arg);
				dim++;
			} else if (arg !== null && typeof arg === 'object') {
				this.parseSettings(arg);
			}
		}, this);

		this.x = this.x || 0;
		this.y = this.y || 0;
	}

	parseSettings(settings) {
		if (settings.dimension && this.dimension < settings.dimension) {
			for (let i = this.dimension; i < settings.dimension; i++) {
				let label = Vector.getDimensionLabel(i);
				this[label] = 0.0;
			}
		}
	}


	/***
		Assign dimension to the n-th dimension in the following order:
		x, y, z, a, b, ... v, w, aa, ab, ...
	***/

	static getDimensionLabel(n) {
		let alphabet = Array(26).fill('').map((n, i) => String.fromCharCode(97 + i)),
			dimSingle = alphabet.slice(-3).concat(alphabet.slice(0, -3));
		return n <= 25 ? dimSingle[n] : alphabet[Math.floor(n / 25) - 1] + alphabet[n % 25];
	}


	/***
		Assign angle to the n-th dimension angle in the greek alphabet:
	***/

	static getAngleLabel(n) {
		let alphabet = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu', 'nu', 'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau', 'upsilon', 'phi', 'chi', 'psi', 'omega'];
		return n <= 23 ? alphabet[n] : alphabet[Math.floor(n / 23) - 1] + alphabet[n % 23];
	}

	


	/******************************************
		Static methods invoked on the class
		e.g. Vector.add(vector1, vector2)
	 ******************************************/

	static sameDimension(vector1, ...vectors) {
		return vectors.filter(vector => vector.dimension !== vector1.dimension).length == 0;
	}

	static add(...vectors) {
		if (!Vector.sameDimension(...vectors)) return;

		let values = Array(vectors[0].dimension).fill(null).map(
			(_, dim) => vectors.reduce((acc, vector) => acc + vector[Vector.getDimensionLabel(dim)], 0)
		);

		return new Vector(...values);
	}

	static addExact(...vectors) {
		if (!Vector.sameDimension(...vectors)) return;

		let values = Array(vectors[0].dimension).fill(null).map(
			(_, dim) => sumExact(vectors.map(vector => vector[Vector.getDimensionLabel(dim)]))
		);

		return new Vector(...values);
	}

	static subtract(vector1, ...vectors) {
		if (!Vector.sameDimension(vector1, ...vectors)) return;

		let values = Array(vector1.dimension).fill(null).map(
			(_, dim) => vectors.reduce((acc, vector) => acc - vector[Vector.getDimensionLabel(dim)], vector1[Vector.getDimensionLabel(dim)])
		);

		return new Vector(...values);
	}

	static subtractExact(...vectors) {
		if (!Vector.sameDimension(...vectors)) return;

		let values = Array(vectors[0].dimension).fill(null).map(
			(_, dim) => sumExact(vectors.map((vector, i) => i == 0 ? vector[Vector.getDimensionLabel(dim)] : - vector[Vector.getDimensionLabel(dim)]))
		);

		return new Vector(...values);
	}

	static mult(vector, scalar) {
		let values = Array(vector.dimension).fill(null).map(
			(_, dim) => vector[Vector.getDimensionLabel(dim)] * scalar
		);

		return new Vector(...values);
	}

	static div(vector, scalar) {
		if (scalar == 0) return;

		let values = Array(vector.dimension).fill(null).map(
			(_, dim) => vector[Vector.getDimensionLabel(dim)] / scalar
		);

		return new Vector(...values);
	}

	static negative(vector) {
		return Vector.mult(vector, -1);
	}

	static dot(vector1, vector2) {
		if (!Vector.sameDimension(vector1, vector2)) return;

		return Array(vector1.dimension).fill(null).map(
			(_, dim) => vector1[Vector.getDimensionLabel(dim)] * vector2[Vector.getDimensionLabel(dim)]
		).reduce((acc, val) => acc + val, 0);
	}

	static length(vector) {
		return Math.sqrt(Vector.dot(vector, vector));
	}

	static normalise(vector) {
		return Vector.div(vector, Vector.length(vector));
	}

	static setLength(vector, scalar) {
		return Vector.mult(Vector.normalise(vector), scalar);
	}

	static limit(vector, scalar) {
		let length = Vector.length(vector);
		return length <= scalar ? vector : Vector.mult(Vector.normalise(vector), scalar);
	}

	static min(vector1, vector2) {
		if (!Vector.sameDimension(vector1, vector2)) return;

		let values = Array(vector.dimension).fill(null).map(
			(_, dim) => Math.min(vector1[Vector.getDimensionLabel(dim)], vector2[Vector.getDimensionLabel(dim)])
		);
		return new Vector(...values);
	}

	static min(vector1, vector2) {
		if (!Vector.sameDimension(vector1, vector2)) return;

		let values = Array(vector.dimension).fill(null).map(
			(_, dim) => Math.max(vector1[Vector.getDimensionLabel(dim)], vector2[Vector.getDimensionLabel(dim)])
		);
		return new Vector(...values);
	}

	static lerp(min, max, fraction) {
		if (!Vector.sameDimension(min, max)) return;

		let newMax = Vector.clone(max);
		return newMax.subtract(min).multiply(fraction).add(min);
	}

	static map(vector, minIn, maxIn, minOut, maxOut) {
		if (!Vector.sameDimension(vector, minIn, maxIn, minOut, maxOut)) return;

		let values = Array(vector.dimension).fill(null).map(
			(_, dim) => {
				let label = Vector.getDimensionLabel(dim);
				return minOut[label] + (vector[label] - minIn[label]) * (maxOut[label] - minOut[label]) / (maxIn[label] - minIn[label]);
			}
		);

		return new Vector(...values);
	}

	static clone(vector) {
		return JSON.parse(JSON.stringify(vector));
	}



	/******************************************
		Public methods invoked on the instance
		which update the instance's coordinates
	 ******************************************/

	get dimension() {
		return Object.keys(this).length;
	}

	update(vector) {
		for (let dim in this) {
			if (vector.hasOwnProperty(dim)) {
				this[dim] = vector[dim];
			}
		}

		return this;
	}

	add(...vectors) {
		vectors.forEach(vector => this.update(Vector.add(this, vector)));
		return this;
	}

	addExact(...vectors) {
		vectors.forEach(vector => this.update(Vector.addExact(this, vector)));
		return this;
	}

	subtract(...vectors) {
		vectors.forEach(vector => this.update(Vector.subtract(this, vector)));
		return this;
	}

	subtractExact(...vectors) {
		vectors.forEach(vector => this.update(Vector.subtractExact(this, vector)));
		return this;
	}

	mult(scalar) {
		return this.update(Vector.mult(this, scalar));
	}

	div(scalar) {
		return this.update(Vector.div(this, scalar));
	}

	negative() {
		return this.update(Vector.negative(this));
	}

	dot(vector) {
		return Vector.dot(this, vector);
	}

	length() {
		return Vector.length(this);
	}

	normalise() {
		return this.update(Vector.normalise(this));
	}

	setLength(scalar) {
		return this.update(Vector.setLength(this, scalar));
	}

	limit(scalar) {
		return this.update(Vector.limit(this, scalar));
	}

	min(vector) {
		return this.update(Vector.min(this, vector));
	}

	min(vector) {
		return this.update(Vector.max(this, vector));
	}

	lerp(min, max, fraction) {
		return this.update(Vector.lerp(min, max, fraction));
	}

	map(minIn, maxIn, minOut, maxOut) {
		return this.update(Vector.map(this, minIn, maxIn, minOut, maxOut));
	}

	clone() {
		return Vector.clone(this);
	}

};

if (typeof exports !== 'undefined') {
	if (typeof module !== 'undefined' && module.exports)
		exports = module.exports = Vector;
	exports.Vector = Vector;
} else {
	window.Vector = Vector;
}