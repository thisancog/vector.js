/*
	parseSettings()
	getDimensionLabel()
	getDimensionFromLabel()
	getAngleLabel()
	clone()
	sameDimension()
	negative()
	cross()
	dot()
	dyadic()
	angle()
	angleToAxis()
	angleX()
	angleY()
	angleZ()
	reduce()
	min()
	max()
	lerp()
	map()
	distance()
	distanceSquared()
	distanceManhattan()
	distanceChebyshev()
	distanceMinkowski()
	scalarTripleProduct()
	vectorTripleProduct()
	matMult()
	.dimension
	.update()

*/


'use strict';

let assert = require('assert'),
	Vector = require('../dist/vector.js');

let vector1 = new Vector(0, 0),
	vector2 = new Vector(12, 21),
	vector3 = new Vector(-5, 9),
	vector4 = new Vector(3, -4),
	vector5 = new Vector(7.4, 1.9),
	vector6 = new Vector(-9.2, 0.01);


// the minimum accuracy, e.g. the maximum difference between the computational and the actual result
let epsilon = 0.00001;


/***************************************************************
	Add and addExact
 ***************************************************************/

describe('add() and addExact()', () => {
	it('should add two 2D vectors regardless of order', () => {
		assert.deepEqual(Vector.add(vector1, vector2), new Vector(12, 21));
		assert.deepEqual(Vector.add(vector2, vector1), new Vector(12, 21));
		assert.deepEqual(Vector.addExact(vector5, vector6), new Vector(-1.8, 1.91));
		assert.deepEqual(Vector.addExact(vector6, vector5), new Vector(-1.8, 1.91));

		// changes vector1!
		vector1.add(vector2);
		assert.deepEqual(vector1, new Vector(12, 21));

		// changes vector5!
		vector5.addExact(vector6);
		assert.deepEqual(vector5, new Vector(-1.8, 1.91));

		vector1 = new Vector(0, 0);
		vector5 = new Vector(7.4, 1.9)
	});

	it('should add several 2D vectors', () => {
		assert.deepEqual(Vector.add(vector1, vector2, vector3), new Vector(7, 30));

		vector2.addExact(vector3, vector5);
		assert.deepEqual(vector2, new Vector(14.4, 31.9));
		vector2 = new Vector(12, 21);
	});

	it('should return input vector if it is the only one', () => {
		assert.deepEqual(Vector.add(vector1), vector1);
	});
});



/***************************************************************
	Subtract and subtractExact
 ***************************************************************/

describe('subtract() and subtractExact()', () => {
	it('should subtract two 2D vectors', () => {
		assert.deepEqual(Vector.subtract(vector1, vector2), new Vector(-12, -21));
		assert.deepEqual(Vector.subtractExact(vector1, vector6), new Vector(9.2, -0.01));
	});

	it('should subtract several 2D vectors', () => {
		assert.deepEqual(Vector.subtract(vector1, vector2, vector3), new Vector(-7, -30));
		assert.deepEqual(Vector.subtractExact(vector1, vector5, vector6), new Vector(1.8, -1.91));

		vector2.subtract(vector3, vector5);
		assert.deepEqual(vector2, new Vector(9.6, 10.1));

		vector5.subtractExact(vector6);
		assert.deepEqual(vector5, new Vector(16.6, 1.89));

		vector2 = new Vector(12, 21);
		vector5 = new Vector(7.4, 1.9);
	});

	it('should return input vector if it is the only one', () => {
		assert.deepEqual(Vector.subtract(vector1), vector1);
		assert.deepEqual(Vector.subtractExact(vector1), vector1);
	});
});



/***************************************************************
	Scalar multiplication and division
 ***************************************************************/

describe('mult() and div()', () => {
	it('should perform scalar multiplication on a 2D vector', () => {
		assert.deepEqual(Vector.mult(vector3, 2), new Vector(-10, 18));

		// changes vector3!
		vector3.mult(2);
		assert.deepEqual(vector3, new Vector(-10, 18));
	});

	it('should perform scalar division on a 2D vector', () => {
		assert.deepEqual(Vector.div(vector3, 2), new Vector(-5, 9));

		// changes vector3!
		vector3.div(2);
		assert.deepEqual(vector3, new Vector(-5, 9));
	});
});



/***************************************************************
	Length, normalise and limit
 ***************************************************************/

describe('length(), normalise(), setLength() and limit()', () => {
	it('should accurately calculate the length/magnitude of a 2D vector', () => {
		assert.deepEqual(Vector.length(vector4), 5);
		assert.deepEqual(vector4.length(), 5);
	});

	it('should calculate the 2D unit vector', () => {
		let result = Vector.normalise(vector2);
		assert.equal(Math.abs(result.x - 0.49613893) < epsilon, true);
		assert.equal(Math.abs(result.y - 0.86824314) < epsilon, true);

		vector2.normalise();
		assert.equal(Math.abs(vector2.x - 0.49613893) < epsilon, true);
		assert.equal(Math.abs(vector2.y - 0.86824314) < epsilon, true);
		vector2 = new Vector(12, 21);
	});

	it("should change the vector's length correctly", () => {
		assert.deepEqual(Vector.limit(vector4, 4), new Vector(2.4, -3.2));
		assert.deepEqual(Vector.setLength(vector4, 5), new Vector(3, -4));

		vector4.limit(4);
		assert.deepEqual(vector4, new Vector(2.4, -3.2));
		vector4.setLength(5);
		assert.deepEqual(vector4, new Vector(3, -4));
	});

	it('should leave the vectors as it is if limit exceeds its length', () => {
		assert.deepEqual(Vector.limit(vector4, 10), new Vector(3, -4));

		vector4.limit(10);
		assert.deepEqual(vector4, new Vector(3, -4));
	});
});




