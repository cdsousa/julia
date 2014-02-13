import Base: getindex, show, +, *
immutable IdentityMatrix{T<:Number} <: AbstractMatrix{T}
	λ::T
end

const I = IdentityMatrix(1)

getindex{T}(I::IdentityMatrix{T}, i::Integer,j::Integer) = ifelse(i==j,one(T),zero(T))

show(io::IO,I::IdentityMatrix) = print(io,"$(typeof(I))\n$(I.λ)*I")

+{T}(B::BitArray{2},I::IdentityMatrix{T}) = +(bitunpack(B),I)
+(I1::IdentityMatrix,I2::IdentityMatrix) = IdentityMatrix(I1.λ+I2.λ)
function +{TA,TI}(A::AbstractMatrix{TA},I::IdentityMatrix{TI})
	n = chksquare(A)
	B = similar(A,promote_type(TA,TI))
	copy!(B,A)
	for i = 1:n
		B[i,i] += I.λ
	end
	B
end
+(I::IdentityMatrix,B::BitArray{2}) = +(IdentityMatrix,bitunpack(B))
+(I::IdentityMatrix,A::AbstractMatrix) = +(A,I)

-(I1::IdentityMatrix,I2::IdentityMatrix) = IdentityMatrix(I1.λ-I2.λ)
-(B::BitArray{2},I::IdentityMatrix) = -(bitunpack(B),I)
-(I::IdentityMatrix,B::BitArray{2}) = -(I,bitunpack(B))
function -{TA,TI<:Number}(A::AbstractMatrix{TA},I::IdentityMatrix{TI})
	n = chksquare(A)
	B = similar(A,promote_type(TA,TI))
	copy!(B,A)
	for i = 1:n
		B[i,i] -= I.λ
	end
	B
end
-(I::IdentityMatrix,B::BitArray{2}) = -(IdentityMatrix,bitunpack(B))
-(I::IdentityMatrix,A::AbstractMatrix) = -(A,I)

*(A::AbstractMatrix,I::IdentityMatrix) = I.λ == 1 ? A : I.λ*A
*(B::BitArray{2},I::IdentityMatrix) = *(bitunpack(B),I::IdentityMatrix)
*(I::IdentityMatrix,B::BitArray{2}) = *(I::IdentityMatrix,bitunpack(B))
*(S::SparseMatrixCSC,I::IdentityMatrix) = I.λ == 1 ? S : I.λ*S
*{Tv,Ti}(I::IdentityMatrix,S::SparseMatrixCSC{Tv,Ti}) = I.λ == 1 ? S : S*I.λ
*(I::IdentityMatrix,A::AbstractMatrix) = I.λ == 1 ? A : I.λ*A
*(A::AbstractMatrix,I::IdentityMatrix) = I.λ == 1 ? A : A*I.λ

*(x::Number,I::IdentityMatrix) = IdentityMatrix(x*I.λ)
