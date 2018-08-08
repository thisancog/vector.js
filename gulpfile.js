var gulp			= require('gulp'),
	compiler		= require('webpack'),
	concat			= require('gulp-concat'),
	minify			= require('gulp-babel-minify'),
	mocha			= require('gulp-mocha'),
	rename			= require('gulp-rename'),
	webpack			= require('webpack-stream');

var scriptSrc		= './src/',
	scriptDst		= './dist/',
	scriptFile		= 'vector.js',
	scriptFileMin	= 'vector.min.js',
	testsSrc = './tests/';



/*****
	Run tests
 *****/

gulp.task('test', function() {
	return gulp.src('./tests/*.js')
			   .pipe(mocha())
});



/*****
	Build and watch
 *****/

gulp.task('build', function() {
	return  gulp.src(scriptSrc + '*.js')
			/*	.pipe(webpack({
						mode: 'development'
					}, compiler))
				.pipe(rename(scriptFile)) */
				.pipe(concat(scriptFile))
				.pipe(gulp.dest(scriptDst))
				.pipe(rename(scriptFileMin))
				.pipe(minify({
					mangle: { keepClassName: true }
				}))
				.pipe(gulp.dest(scriptDst));
});

gulp.task('default', gulp.series(['build']), function() {
	gulp.watch(scriptSrc + '*.js', ['build']);
});


