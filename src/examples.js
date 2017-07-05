
const simple =
`
console.log('BlueMath version',bluemath.version);
let A = new bluemath.NDArray([
  [2,0],
  [5,1]
]);
console.log('A',A.toString());
`;

const identity_and_multiplication =
`
const bm = bluemath;

let A = bm.eye(3,'i32');
console.log('Identity',A);
let B = bm.mul(A,5);
let C = bm.mul(A,2);
console.log('B',B);
console.log('C',C);
let D = bm.mul(B,C);
console.log('B*C',D);
`;

const determinant = 
`
const {NDArray,linalg} = bluemath;

let A = new NDArray([
    [3,4,5],
    [0,3,4],
    [1,3,5]
]);

console.log('Determinant = ',linalg.det(A));
`;

let BMSHELL_EXAMPLES = {
  'Simple' : simple,
  'Identity and Multiplication' : identity_and_multiplication,
  'Determinant' : determinant
};