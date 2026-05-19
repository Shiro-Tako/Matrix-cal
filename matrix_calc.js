const EPS = 1e-8;

function clone(A){ return A.map(r=>[...r]); }
function eye(n){ return Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i===j?1:0)); }
function mul(A,B){
  const n=A.length,m=B[0].length,k=B.length;
  const C=Array.from({length:n},()=>Array(m).fill(0));
  for(let i=0;i<n;i++) for(let j=0;j<m;j++) for(let t=0;t<k;t++) C[i][j]+=A[i][t]*B[t][j];
  return C;
}
function transpose(A){ return A[0].map((_,j)=>A.map(r=>r[j])); }
function norm(v){ return Math.sqrt(v.reduce((s,x)=>s+x*x,0)); }

function qrDecompose(A){
  const n=A.length, Q=Array.from({length:n},()=>Array(n).fill(0)), R=Array.from({length:n},()=>Array(n).fill(0));
  const V=transpose(A), U=[];
  for(let i=0;i<n;i++){
    let u=[...V[i]];
    for(let j=0;j<i;j++){
      const dot=V[i].reduce((s,x,k)=>s+x*U[j][k],0);
      R[j][i]=dot;
      u=u.map((x,k)=>x-dot*U[j][k]);
    }
    const un=norm(u)||EPS;
    R[i][i]=un;
    U.push(u.map(x=>x/un));
  }
  for(let i=0;i<n;i++) for(let j=0;j<n;j++) Q[j][i]=U[i][j];
  return {Q,R};
}

function qrEigenvalues(A, it=220){
  let Ak=clone(A);
  for(let i=0;i<it;i++){
    const {Q,R}=qrDecompose(Ak);
    Ak=mul(R,Q);
  }
  return Ak.map((r,i)=>r[i]);
}

function rref(M){
  const A=clone(M); const rows=A.length, cols=A[0].length;
  let r=0; const pivots=[];
  for(let c=0;c<cols && r<rows;c++){
    let p=r;
    for(let i=r+1;i<rows;i++) if(Math.abs(A[i][c])>Math.abs(A[p][c])) p=i;
    if(Math.abs(A[p][c])<EPS) continue;
    [A[r],A[p]]=[A[p],A[r]];
    const div=A[r][c]; for(let j=c;j<cols;j++) A[r][j]/=div;
    for(let i=0;i<rows;i++) if(i!==r){
      const f=A[i][c]; if(Math.abs(f)<EPS) continue;
      for(let j=c;j<cols;j++) A[i][j]-=f*A[r][j];
    }
    pivots.push(c); r++;
  }
  return {A,pivots};
}

function eigenvectorsForLambda(A,l){
  const n=A.length;
  const M=A.map((row,i)=>row.map((v,j)=>v-(i===j?l:0)));
  const {A:R,pivots}=rref(M);
  const free=[]; for(let c=0;c<n;c++) if(!pivots.includes(c)) free.push(c);
  const vecs=[];
  for(const f of free){
    const x=Array(n).fill(0); x[f]=1;
    for(let i=pivots.length-1;i>=0;i--){
      const pc=pivots[i]; let s=0;
      for(let j=pc+1;j<n;j++) s+=R[i][j]*x[j];
      x[pc]=-s;
    }
    const xn=norm(x)||1;
    vecs.push(x.map(v=>v/xn));
  }
  return vecs;
}

function inverse(A){
  const n=A.length; const M=A.map((row,i)=>[...row,...eye(n)[i]]);
  const {A:R,pivots}=rref(M);
  if(pivots.length<n) return null;
  return R.map(row=>row.slice(n));
}

function uniqueApprox(vals){
  const out=[];
  for(const v of vals) if(!out.some(x=>Math.abs(x-v)<1e-5)) out.push(v);
  return out;
}

function fmtMat(A){ return A.map(r=>`[ ${r.map(v=>v.toFixed(6)).join(' , ')} ]`).join('\n'); }

export function buildReport(A){
  const n=A.length;
  if(![2,3].includes(n) || A.some(r=>r.length!==n)) throw new Error('Matrix must be 2x2 or 3x3');
  const vals = qrEigenvalues(A);
  const uniq = uniqueApprox(vals);

  let report = `A =\n${fmtMat(A)}\n\n1) Eigenvalues (approx)\n`;
  report += vals.map((v,i)=>`λ${i+1} = ${v.toFixed(6)}`).join('\n') + '\n\n2) Eigenvectors\n';

  const eigenPairs=[];
  for(const l of uniq){
    const vecs=eigenvectorsForLambda(A,l);
    report += `For λ = ${l.toFixed(6)}\n`;
    if(!vecs.length) report += '  - Could not compute a real eigenvector numerically.\n';
    vecs.forEach((v,idx)=>{ report += `  v${idx+1} = [${v.map(x=>x.toFixed(6)).join(', ')}]\n`; eigenPairs.push({l,v}); });
  }

  report += '\n3) Diagonalizable?\n';
  if(eigenPairs.length < n){
    report += 'No\n\n4) Cannot construct P and P^-1\n5) Reason: number of independent eigenvectors is less than matrix size.\n';
    return report;
  }

  const P = Array.from({length:n},(_,r)=>eigenPairs.slice(0,n).map(p=>p.v[r]));
  const PInv = inverse(P);
  if(!PInv){
    report += 'No\n\n4) Cannot construct P and P^-1\n5) Reason: eigenvector matrix P is singular.\n';
    return report;
  }

  const D = mul(PInv,mul(A,P));
  report += 'Yes\n\n4) P, P^-1, D = P^-1 A P\n';
  report += `P =\n${fmtMat(P)}\n\nP^-1 =\n${fmtMat(PInv)}\n\nD =\n${fmtMat(D)}\n\n5) If not diagonalizable, reason is shown above.\n`;
  return report;
}
