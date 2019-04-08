

#ifndef BDTS_MATH__INCLUDED
#define BDTS_MATH__INCLUDED

#include "bdts_base.h"
#include <math.h>


#if BDTS_IS_CPP_COMPILE 
namespace bdts
{
#endif
	const float bdts_pi = 3.14159f;

	static float bdts_max(float a, float b)
	{
		return a < b ? b : a;
	}

	static float bdts_min(float a, float b)
	{
		return a < b ? a : b;
	}


	// log natural aka 1/x integral
	static float logN(float x)
	{
		return logf(x) / logf(2.7182f);
	}


	inline float rndFromIndex(unsigned long index, unsigned long offset = 0, unsigned long optSeed = 0)
	{

		unsigned long stuff = (index + 1) * 179443609 + (offset + 1) * 198510517 + (optSeed + 1) * 15485863;
		stuff = stuff % 1628149;
		static const unsigned long floatConv = 7919;
		stuff = stuff % floatConv;
		return   (float)(stuff / ((double)floatConv));
	}




	static float tosign(float x) {	return x < 0.0f ? -1.0f : 1.0f;	}

		struct float2
		{
			inline float2() { x = 0.0f; y = 0.0f; }
			inline explicit float2(float xx) { x = xx; y = xx; }
			inline float2(float xx, float yy) { x = xx; y = yy; }
			inline float2(const float2& v) { this->x = v.x; this->y = v.y; }
			float2& operator=(float2 v)
			{
				this->x = v.x;
				this->y = v.y;
				return *this;
			}

			inline float2 GetSign() const { return float2(tosign(x), tosign(y)); }
			inline float2 GetCross() const { return float2(-y, x); }
			inline float GetAngle() const { return atan2f(y, x); }
			inline float2 operator- () const { return float2(-x,-y); } // unary
			inline float2 operator+ (float2 v2) const { return float2(x + v2.x, y + v2.y); }
			inline float2 operator- (float2 v2) const { return float2(x - v2.x, y - v2.y); }
			inline float2 operator* (float2 v2) const { return float2(x * v2.x, y * v2.y);}
			inline float2 operator/ (float2 v2) const { return float2(x / v2.x, y / v2.y); }
			inline float2 operator* (float  v2) const { return float2(x * v2, y * v2);}
			inline float2 operator/ (float  v2) const { return float2(x / v2, y / v2); }

			inline float2& operator += (float2 v) { *this = *this + v; return *this; }
			inline float2& operator -= (float2 v) { *this = *this - v; return *this; }
			inline float2& operator *= (float2 v) { *this = *this * v; return *this; }
			inline float2& operator /= (float2 v) { *this = *this / v; return *this; }

			inline friend float2 operator* (float   v1, float2 v2) { return float2(v1) * v2; }
			inline friend float2 operator/ (float   v1, float2 v2) { return float2(v1) / v2; }
			float x;
			float y;
		};

		__declspec(align(16)) struct float3
		{
			inline float3() { x = 0.0f; y = 0.0f; z = 0.0f; }
			inline explicit float3(float xx) { x = xx; y = xx; z = xx; }
			inline float3(float xx, float yy, float zz) { x = xx; y = yy; z = zz; }
			inline float3(const float3& v) { this->x = v.x; this->y = v.y; this->z = v.z; }
			float3& operator=(float3 v)
			{
				this->x = v.x;
				this->y = v.y;
				this->z = v.z;
				return *this;
			}

			inline float3 GetSign() const { return float3(tosign(x), tosign(y), tosign(z)); }
			
			inline float3 operator- () const { return float3(-x, -y, -z); } // unary
			inline float3 operator+ (float3 v2) const { return float3(x + v2.x, y + v2.y, z + v2.z); }
			inline float3 operator- (float3 v2) const { return float3(x - v2.x, y - v2.y, z - v2.z);}
			inline float3 operator* (float3 v2) const { return float3(x * v2.x, y * v2.y, z * v2.z);}
			inline float3 operator/ (float3 v2) const { return float3(x / v2.x, y / v2.y, z / v2.z); }
			inline float3 operator* (float  v2) const { return *this * v2; }
			inline float3 operator/ (float  v2) const { return *this / v2; }

			inline float3& operator += (float3 v) { *this = *this + v; return *this; }
			inline float3& operator -= (float3 v) { *this = *this - v; return *this; }
			inline float3& operator *= (float3 v) { *this = *this * v; return *this; }
			inline float3& operator /= (float3 v) { *this = *this / v; return *this; }

			inline friend float3 operator* (float   v1, float3 v2) { return float3(v1) * v2; }
			inline friend float3 operator/ (float   v1, float3 v2) { return float3(v1) / v2; }
			float x;
			float y;
			float z;
			float w;// unused
		};


		__declspec(align(16)) struct float4
		{
			inline float4() { x = 0.0f; y = 0.0f; z = 0.0f; w = 0.0f; }
			inline explicit float4(float xx) { x = xx; y = xx; z = xx; w = xx; }
			inline float4(float xx, float yy, float zz, float ww) { x = xx; y = yy; z = zz; w = ww; }
			inline float4(const float4& v) {
				this->x = v.x;
				this->y = v.y;
				this->z = v.z;
				this->w = v.w;
			}

			float4& operator=(float4 v)
			{
				this->x = v.x;
				this->y = v.y;
				this->z = v.z;
				this->w = v.w;
				return *this;
			}

			inline float4 GetSign() const { return float4(tosign(x), tosign(y), tosign(z),tosign(w)); }

			inline float4 operator- () const { return float4(-x, -y, -z,-w); } // unary
			inline float4 operator+ (float4 v2) const { return float4(x + v2.x, y + v2.y, z + v2.z, w + v2.w); }
			inline float4 operator- (float4 v2) const { return float4(x - v2.x, y - v2.y, z - v2.z, w - v2.w); }
			inline float4 operator* (float4 v2) const { return float4(x * v2.x, y * v2.y, z * v2.z, w * v2.w); }
			inline float4 operator/ (float4 v2) const { return float4(x / v2.x, y / v2.y, z / v2.z, w / v2.w); }
			inline float4 operator* (float  v2) const { return float4(x * v2, y * v2, z* v2, w * v2);}
			inline float4 operator/ (float  v2) const { return float4(x /v2, y / v2, z /v2, w / v2);}

			inline float4& operator += (float4 v) { *this = *this + v; return *this; }
			inline float4& operator -= (float4 v) { *this = *this - v; return *this; }
			inline float4& operator *= (float4 v) { *this = *this * v; return *this; }
			inline float4& operator /= (float4 v) { *this = *this / v; return *this; }

			inline friend float4 operator* (float   v1, float4 v2) { return float4(v1) * v2; }
			inline friend float4 operator/ (float   v1, float4 v2) { return float4(v1) / v2; }
			float x;
			float y;
			float z;
			float w;
		};


		__declspec(align(16)) class Matrix22
		{
		public:
			inline Matrix22() { m_mat = float4(); }
			inline explicit Matrix22(const float4& xyzw) { m_mat = xyzw; }
			inline Matrix22(const Matrix22& m) { m_mat = m.m_mat; }
			//inline explicit Matrix22(EIdentityTag) { m_mat = Vector4(1.0f, 0.0f, 0.0f, 1.0f); }
			//inline explicit Matrix22(EZeroTag) { m_mat = Vector4(kZero); }

			inline void Set(const float4& xyzw) { m_mat = xyzw; }
			inline float4 Get() const { return m_mat; }

			inline float2 operator* (const float2& vec) const {
				return  float2(m_mat.x, m_mat.z) * vec.x + float2(m_mat.y, m_mat.w) * vec.y;
			}

			inline Matrix22 operator* (const Matrix22& mat) const
			{
				float2 col0 = (*this) * float2(mat.Get().x, mat.Get().z);
				float2 col1 = (*this) * float2(mat.Get().y, mat.Get().w);
				return Matrix22(float4(col0.x, col1.x,
					col0.y, col1.y));
			}

			inline Matrix22 operator* (const float mult) const
			{
				return Matrix22(m_mat * mult);
			}

			static Matrix22 rotation(float rot)
			{
				return Matrix22({ cosf(rot) , -sinf(rot), sinf(rot) , cosf(rot) });
			}

			inline Matrix22 operator+ (Matrix22 v2) const { return Matrix22(m_mat + v2.m_mat); }
			inline Matrix22 operator- (Matrix22 v2) const { return Matrix22(m_mat - v2.m_mat); }

		private:

			// x y
			// z w
			float4 m_mat;
		};




		struct Line2d
		{
			Line2d(float2 inpi, float2 inpf) :
				pi(inpi), pf(inpf)
			{
			}
			float2 pi;// init
			float2 pf;// final
			bool hasLength() const;
			float2 getVec() const;
			// y = a x + b
			float getCompA() { float2 diff = pf - pi; return diff.y / diff.x; }
			float getCompB() { return  pi.y - (getCompA() * pi.x); }


			bool isNearlyVertical() { return fabsf(getCompA()) > 1000.0f; }
			bool isNearlyHorizontal() { return fabsf(getCompA()) < 0.001f; }
			float fromXToY(float x, float failValueY)
			{
				if (isNearlyVertical())
				{
					return failValueY;
				}

				float a = getCompA();
				float b = getCompB();
				return a * x + b;
			}

			float fromYToX(float y, float failValueX)
			{
				if (isNearlyVertical())
				{
					return failValueX;
				}

				if (isNearlyHorizontal())
				{
					return failValueX;
				}

				float a = getCompA();
				float b = getCompB();
				return (y - b) / a;
			}

			// line point distance mins
			bool getStartStopForDistance(float distance, float2& outStart, float2& outStop, bool limited = true);


			float dirXToCoordX(float dirX);

			// the pi to pf forms a vector which spits the plane. The cross is the up direction
			bool isUpRelative(const float2& ptest) const;
			float distToLine(const float2& ptest) const;
			bool isPointAlongLine(const float2& ptest) const;
			float projPointLineFloat(const float2& ptest) const;
			float length() const;
		};





		struct Triangle2d
		{
			Triangle2d()
			{
				v[0] = float2();
				v[1] = float2();
				v[2] = float2();
				vClip[0] = true;
				vClip[1] = true;
				vClip[2] = true;
			}

			Triangle2d(float2 v0, float2 v1, float2 v2)
			{
				v[0] = v0;
				v[1] = v1;
				v[2] = v2;
				vClip[0] = true;
				vClip[1] = true;
				vClip[2] = true;
			}
			float2 v[3];
			bool vClip[3];
			bool isLeft() const;
			bool hasVolume() const;
			float getArea() const;
			float robustness() const;
			Triangle2d toRight() const;
			Triangle2d toLeft() const;
			float2 getCentroid() const;
			float getMaxAngle() const;
		};
		static float2 horzNormDir(float2 in)
		{
			if (in.y <= 0.0f)
			{
				in.x = (in.x < 0.0f ? -0.999f : 0.999f);
			}
			return in;
		}

#if BDTS_IS_CPP_COMPILE 
}
#endif

#endif //BDTS_MATH__INCLUDED