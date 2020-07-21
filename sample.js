const Fresnel = require("./index.js");

const L = 200;

console.log('<svg viewBox="0 0 200 200" xmlns="http://www.w3.org/2000/svg">');
console.log('<polyline fill="none" stroke="black" points="');
for (let i = 0; i < 200; i++) {
  var pt = Fresnel(i / 60);
  console.log((pt.x*L + L/10) + "," + (pt.y*L + L/10))
}
console.log('"/>');
console.log('</svg>')
