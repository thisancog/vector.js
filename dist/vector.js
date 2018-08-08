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


"use strict";

class Vector2D {
	constructor(x, y) {
		this.x = parseFloat(x) || 0;
		this.y = parseFloat(y) || 0;
	}


	/******************************************
		Static methods invoked on the class
		e.g. Vector2D.add(vector1, vector2)
	 ******************************************/

	static add(...vectors) {
		let x = vectors.reduce((acc, vector) => acc + vector.x, 0),
			y = vectors.reduce((acc, vector) => acc + vector.y, 0);
		return new Vector2D(x, y);
	}

	static addExact(...vectors) {
		let x = sumExact(vectors.map((vector) => vector.x)),
			y = sumExact(vectors.map((vector) => vector.y));
		return new Vector2D(x, y);
	}

	static subtract(vector1, ...vectors) {
		let x = vectors.reduce((acc, vector) => acc - vector.x, vector1.x),
			y = vectors.reduce((acc, vector) => acc - vector.y, vector1.y);
		return new Vector2D(x, y);
	}

	static subtractExact(...vectors) {
		let x = sumExact(vectors.map((vector, i) => i === 0 ? vector.x : -vector.x)),
			y = sumExact(vectors.map((vector, i) => i === 0 ? vector.y : -vector.y));
		return new Vector2D(x, y);
	}

	static mult(vector, scalar) {
		return new Vector2D(vector.x * scalar, vector.y * scalar);
	}

	static div(vector, scalar) {
		return new Vector2D(vector.x / scalar, vector.y / scalar);
	}

	static length(vector) {
		return Math.sqrt(vector.x * vector.x + vector.y * vector.y);
	}

	static normalise(vector) {
		return Vector.div(vector, Vector2D.length(vector));
	}

	static setLength(vector, scalar) {
		return Vector.mult(Vector2D.normalise(vector), scalar);
	}

	static limit(vector, scalar) {
		let length = Vector2D.length(vector);
		return length <= scalar ? vector : Vector2D.mult(Vector2D.normalise(vector), scalar);
	}



	/******************************************
		Public methods invoked on the instance
		which update the instance's coordinates
	 ******************************************/

	update(vector) {
		this.x = vector.x;
		this.y = vector.y;
	}

	add(...vectors) {
		vectors.forEach((vector) => this.update(Vector2D.add(this, vector)));
	}

	addExact(...vectors) {
		vectors.forEach((vector) => this.update(Vector2D.addExact(this, vector)));
	}

	subtract(...vectors) {
		vectors.forEach((vector) => this.update(Vector2D.subtract(this, vector)));
	}

	subtractExact(...vectors) {
		vectors.forEach((vector) => this.update(Vector2D.subtractExact(this, vector)));
	}

	mult(scalar) {
		this.update(Vector2D.mult(this, scalar));
	}

	div(scalar) {
		this.update(Vector2D.div(this, scalar));
	}

	length() {
		return Vector2D.length(this);
	}

	normalise() {
		this.update(Vector2D.normalise(this));
	}

	setLength(scalar) {
		this.update(Vector2D.setLength(this, scalar));
	}

	limit(scalar) {
		this.update(Vector2D.limit(this, scalar));
	}

};

if (typeof exports !== 'undefined') {
	if (typeof module !== 'undefined' && module.exports)
		exports = module.exports = Vector2D;
	exports.Vector2D = Vector2D;
} else {
	window.Vector2D = Vector2D;
}

