#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <vector>

#ifdef __cplusplus  
extern "C" {  
#endif
#include <pbm.h>
#ifdef __cplusplus
}
#endif

#define CURVE_BUF_SIZE 200

struct LL {
    char degree;
    char minute;
    float second;
    void toString() {
        printf("%3d°%02d'%02f\"", degree, minute, second);
    }
};

struct LatitudeLongitudeSpeed {
    struct LL latitude;
    struct LL longitude;
    float speed;
};

static struct LatitudeLongitudeSpeed * data = 0;
static int lines = 0;

/**
 * Graph will draw the a 2d function to 2d PBM.
**/
template <typename T>
class Graph {
private:
    T mMax;
    T mMin;
    T mThreshold;
    int mSize;
    float mScaleY;
    float mScaleX;
    T * mData;
    bit ** mBits;
    int mWidth;
    int mHeight;
    int * mVerticalLines;

public:
    Graph() : mMax(0), mMin(0), mThreshold(0), mSize(0), mScaleX(1), mScaleY(1), mData(NULL),
        mBits(NULL), mWidth(0), mHeight(0), mVerticalLines(NULL) {
    }

private:
    void findMaxDiff() {
        T max, min;
        max = min = mData[0];
        for (int i = 1; i < mSize; i++) {
            T d = mData[i];
            if (d > max)
                max = d;
            else if (d < min)
                min = d;
        }
        mMax = max;
        mMin = min;
        if (mMin > 0)
            mMin = 0;

        if (mMax < 0)
            mMax = 0;
        //printf("mMax = %f, mMin = %f\n", mMax, mMin);
    }

public:
    void setMax(int max) {
        mMax = max;
    }

    void setThreshold(T threshold) {
        mThreshold = threshold;
    }

    void setVerticalLines(int * lines) {
        mVerticalLines = lines;
    }

    void setData(T* data, int size) {
        mMax = 0;
        mMin = 0;
        mData = data;
        mSize = size;
        if (data == NULL)
            return;
        findMaxDiff();
    }

private:
    inline void validate(int& x, int& y) {
        if (x >= mWidth)
            x = mWidth - 1;
        if (y >= mHeight)
            y = mHeight - 1;
        if (x < 0)
            x = 0;
        if (y < 0)
            y = 0;
    }

    void drawLine(int x0, int y0, int x1, int y1) {
        if (mBits == NULL)
            return;

        int dx = x1 - x0;
        int dy = y1 - y0;
        if (dx == 0 && dy == 0)
            return;

        int startX, endX;
        int startY, endY;
        if (x0 < x1) {
            startX = x0;
            endX = x1;
        } else {
            startX = x1;
            endX = x0;
        }
        if (y0 < y1) {
            startY = y0;
            endY = y1;
        } else {
            startY = y1;
            endY = y0;
        }

        validate(startX, startY);
        validate(endX, endY);

        // Virtical
        if (dx == 0) {
            int x = x0;
            int y = y0;
            validate(x, y);
            for (int j = startY; j <= endY; j++) {
                mBits[j][x] = 1;
            }
            return;
        }

        // Horizontal
        if (dy == 0) {
            int x = x0;
            int y = y0;
            validate(x, y);
            for (int i = startX; i <= endX; i++) {
                mBits[y][i] = 1;
            }
            return;
        }

        // Draw in two direction for avoiding drawing a discontinue line.

        // y = m(x-x0) + y0
        float slopForY = (float)dy / (float)dx;
        for (int i = startX; i <= endX; i++) {
            int j = static_cast<int>(slopForY * (i - x0) + y0);
            if (j < 0 || j >= mHeight)
                continue;
            mBits[j][i] = 1;
        }

        // x = m(y-y0) + x0
        float slopForX = (float)dx / (float)dy;
        for (int j = startY; j <= endY; j++) {
            int i =  static_cast<int>(slopForX * (j - y0) + x0);
            if (i < 0 || i >= mWidth)
                continue;
            mBits[j][i] = 1;
        }
        return;
    }

    void drawMarkedLines() {
        int center = mMax * mScaleY;
        drawLine(0, center, mWidth, center);
        if (mThreshold != 0) {
            int t =  mThreshold * mScaleY;
            drawLine(0, center + t, mWidth, center + t);
            drawLine(0, center - t, mWidth, center - t);
        }
        if (mVerticalLines != NULL) {
            for (int i = 0; i < CURVE_BUF_SIZE; i++) {
                if (mVerticalLines[i] <= 0)
                    return;
                int x = 0;
                if (mScaleX > 0)
                    x = static_cast<int>((mVerticalLines[i] + 0.5f) * mScaleX);
                if (mScaleX == 0)
                    x = mVerticalLines[i];
                
                drawLine(x, 0, x, mHeight - 1);
            }
        }
    }

    void dataToPBM() {
        int center = mMax * mScaleY;
        if (mData == NULL)
            return;
        if (mScaleX < 1) {
            printf("Not support the mScaleX < 1 now\n");
        } if (mScaleX == 1) {
            //printf("mScaleX = 1\n");
            drawMarkedLines();
            for (int i = 0; i < mSize; i++) {
                int x = i;
                int y = (mMax - mData[i]) * mScaleY;
                drawLine(x, y, x, center);
            }
        } else {
            // mScaleX > 1
            //printf("mScaleX > 1\n");
            drawMarkedLines();
            int lastX = 0, lastY = center;
            for (int i = 0; i < mSize; i++) {
                // draw in center of the gap (mScale)
                int x = static_cast<int>((i + 0.5f) * mScaleX);
                int y = (mMax - mData[i]) * mScaleY;
                drawLine(lastX, lastY, x, y);
                lastX = x;
                lastY = y;
            }
        }
    }

public:
    void write(const char * name, int w, int h) {
        mWidth = w;
        mHeight = h;
        mScaleY = h / (mMax - mMin);
        mScaleX = w / mSize;
        bit ** bits = pbm_allocarray(w, h);
        for (int i = 0; i < w; i++)
            for (int j = 0; j < h; j++)
                bits[j][i] = 0;
        mBits = bits;

        FILE * file = pm_openw(name);
        if (file == NULL)
            return;
        dataToPBM();
        pbm_writepbm(file, bits, w, h, true);
        printf("Graph File:%s Max: %f Min: %f\n", name, (float)mMax, (float)mMin);

        pm_close(file);
        pbm_freearray(mBits, h);
        mBits = NULL;
    }
};


class Gaussian {
public:
    const double * mData;
    double * mResult;
    int mFilterSize;
    const int mSize;
    Gaussian(double * data, size_t s) : mData(data), mResult(NULL), mSize(s) {
    }

    ~Gaussian() {
        if (mResult == NULL)
            delete [] mResult;
    }

    void filter(size_t filterSize) {
        if (mResult == NULL)
            mResult = new double [mSize];
        mFilterSize = filterSize;
        int half = filterSize / 2;
        for (int i = 0; i < half; i++) {
            mResult[i] = mData[i];
        }
        for (int i = mSize - half; i < mSize; i++) {
            mResult[i] = mData[i];
        }
        for (int i = half; i < mSize - half; i++) {
            double sum = 0;
            for (int f = -half; f < (int)filterSize - half; f++) {
                sum += mData[f + i];
            }
            double avg = sum / (int)filterSize;
            mResult[i] = avg;
        }
    }

    void writeNameByFilter() {
        if (mResult == NULL)
            return;
        char name[30];
        sprintf(name, "/tmp/pbm/gau_%d.pbm", mFilterSize);
        Graph<double> g;
        g.setData(mResult, mSize);
        g.write(name, 1200, 300);
    }

    void write(const char * name) {
        if (mResult == NULL)
            return;
        Graph<double> g;
        g.setData(mResult, mSize);
        g.write(name, 1200, 300);
    }
};

class Wavelet {
public:
    double * mData;
    float * mCenter;
    unsigned int * mShift;
    int mSize;
    size_t mRawSize;
    unsigned int mLevel;

    Wavelet(size_t s) : mData(NULL), mSize(0), mRawSize(s) {
        int exp = 0;
        do {
            exp++;
            s >>= 1;
        } while (s > 1);
        mLevel = exp;
        mSize = pow(2, mLevel + 1);
        mData = new double [mSize];
        mShift = new unsigned int [mLevel];
        mCenter = new float [mSize];
        printf ("wavelet Level=%d, Size=%u\n", mLevel, mSize);

        for (unsigned int i = 0; i < mLevel; i++) {
            mShift[i] = 0;
        }
    }

    ~Wavelet() {
        delete [] mData;
        delete [] mShift;
    }

    void setShift(unsigned int level, unsigned int shift) {
        if (shift != 0 && shift != 1) {
            printf("shift can only be 0 or 1.\n");
            return;
        }
        mShift[level] = shift;
    }

    void setData(double * data, unsigned int hideLevel, unsigned int level = 0) {
#define SHOW_TEXT false
        double * resMinus = mData, * res_head = mData;
        double * ptr = data;
        double * resPlus_head = new double [mSize / 2];
        double * resPlus = resPlus_head;
        if (level == 0)
            level = mLevel;
        for (unsigned int l = 0, size = mRawSize;
                size != 0 && l < level;
                l++) {
            if (l >= hideLevel && SHOW_TEXT)
                printf("--- level %2d ---\n", l);
            for (unsigned int j = mShift[l]; j < size - 1; j += 2, resMinus++, resPlus++) {
                *resMinus = (ptr[j + 1] - ptr[j]);
                *resPlus = (ptr[j + 1] + ptr[j]);
                if (l >= hideLevel && SHOW_TEXT)
                    printf("%d\t%F\t%F\n", j, *resMinus, *resPlus);
            }
            if (l >= hideLevel) {
                int datasize = (size - 1 - mShift[l]) / 2;
                const int fileN = 30;
                char plusname[fileN], minusname[fileN];
                snprintf(plusname, fileN - 1, "/tmp/pbm/plus_%02d.pbm", l);
                snprintf(minusname, fileN - 1, "/tmp/pbm/minus_%02d.pbm", l);
                Graph<double> plus, minus;
                plus.setData(resPlus_head, datasize);
                plus.write(plusname, 600, 300);

                minus.setData(res_head, datasize);
                minus.write(minusname, 600, 300);

                if (SHOW_TEXT)
                    printf("\n");
            }
            ptr = resPlus_head;
            resPlus = resPlus_head;
            res_head += size / 2;
            resMinus = res_head;
            size = size / 2;
        }
        delete [] resPlus_head;
    }
};


static void toDegree(char * ll, struct LL& llstruct) {
    if (ll == NULL)
        return;
    char buf[10];
    char * ptr = ll, * idx = strrchr(ll, '.');
    if (idx == NULL)
        return;
    int degreeSize = idx - ptr - 2;
    for (int i = 0; i < degreeSize; i++)
        buf[i] = ptr[i];
    buf[degreeSize] = 0;
    llstruct.degree = atoi(buf);

    ptr += degreeSize;
    buf[0] = ptr[0];
    buf[1] = ptr[1];
    buf[2] = 0;
    llstruct.minute = atoi(buf);

    ptr += 2;
    llstruct.second = atof(ptr) * 60;

    //printf("%s %d %d %f\n", ll, llstruct.degree, llstruct.minute, llstruct.second);
}

static void readFile(FILE* fp) {
    if (fp == NULL)
        return;
    data = new struct LatitudeLongitudeSpeed[2000];
    if (data == NULL)
        return;

    struct LatitudeLongitudeSpeed * ptr = data;
    do {
        char latbuf[15] = {0}, lonbuf[15] = {0};
        int ret = fscanf(fp, "%s\t%s\t%f\n", latbuf, lonbuf, &ptr->speed);
        if (ret == EOF) {
            break;
        }
        toDegree(latbuf, ptr->latitude);
        toDegree(lonbuf, ptr->longitude);
        lines++;
    } while (++ptr < data+2000);
}

static float Q_rsqrt(float number) {
	long i;
	float x2, y;
	const float threehalfs = 1.5F;
 
	x2 = number * 0.5F;
	y  = number;
	i  = * ( long * ) &y;                       // evil floating point bit level hacking
	i  = 0x5f3759df - ( i >> 1 );
	y  = * ( float * ) &i;
	y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
//      y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed
 
	return y;
}

class Vector {
public:
    float mX;
    float mY;
    float mNorm;

    Vector() : mX(0), mY(0), mNorm(0) {
    }

    Vector(float x, float y) : mX(x), mY(y), mNorm(Vector::norm()) {
    }

    Vector(const Vector& v) : mX(v.mX), mY(v.mY), mNorm(v.mNorm) {
    }

    void set(float x, float y) {
        mX = x;
        mY = y;
        mNorm = Vector::norm();
    }

private:
    inline float norm() {
        //return 1 / Q_rsqrt(mX*mX + mY*mY);
        return sqrt(((double)mX)*mX + ((double)mY)*mY);
    }

public:
    inline static float innerProduct(const Vector& a, const Vector& b) {
        return a.mX * b.mX + a.mY * b.mY;
    }

    inline static float crossProduct(const Vector& a, const Vector& b) {
        return a.mX * b.mY - b.mX * a.mY;
    }

    inline static float cos0(const Vector& a, const Vector& b) {
        return innerProduct(a, b) / a.mNorm / b.mNorm;
    }

    inline static float sin0(const Vector& a, const Vector& b) {
        return crossProduct(a, b) / a.mNorm / b.mNorm;
    }

    inline static float distance(const Vector& a, const Vector& b) {
        float dx = a.mX - b.mX;
        float dy = a.mY - b.mY;
        return sqrt(dx * dx + dy * dy);
    }

    void toString() {
        printf("%f,%f,%f", mX, mY, mNorm);
    }
};

class FindCurve {
private:
    double * mData;
    int * mCurveLocation;
    ssize_t mSize;
    double mThreshold;
    double mThreshold2;

public:
    class Curve {
    public:
        double * mData;
        size_t mSize;
        double mTotal;
        unsigned int mPeak;
        Curve() : mData(NULL), mSize(0), mTotal(0), mPeak(0) {
        }

        Curve(Curve & o) : mData(o.mData), mSize(o.mSize) {
            o.mData = NULL;
            o.mSize = 0;
        }

        ~Curve() {
            if (mData != NULL)
                delete [] mData;
        }

        bool setData(const double * data, size_t size) {
            if (size < 4)
                return false;
            mSize = size;
            mData = new double [mSize];
            int sign = data[0] > 0 ? 1 : -1;
            mTotal = 0;
            for (int i = 0; i < size; i++) {
                double angle = data[i] * sign;
                mTotal += angle;
                mData[i] = angle;
            }
            if (mTotal < 20)
                return false;
            double peak = mTotal / 2.0;
            double sum = 0;
            for (int i = 0; i < size; i++) {
                sum += mData[i];
                if (sum >= peak) {
                    mPeak = i;
                    break;
                }
            }
            return true;
        }

        void writeSerial(unsigned int s) {
            if (mData == NULL)
                return;
            char name[30];
            snprintf(name, 30, "/tmp/pbm/curve_%04u.pbm", s);
            Graph<double> g;
            g.setData(mData, mSize);
            int vlines[2] = {(int)mPeak, 0};
            g.setVerticalLines(vlines);
            g.write(name, 100, 300);
            printf("Size %lu, Total %6f Peak at %d\n", mSize, mTotal, mPeak);
        }

        double operator[](int i) {
            if (mData == NULL)
                return 0;
            return mData[i];
        }
    };

    FindCurve(double * data, ssize_t size) :
        mData(data), mCurveLocation(NULL), mSize(size), mThreshold(2), mThreshold2(1) {
        mCurveLocation = new int [CURVE_BUF_SIZE];
        for (int i = 0; i < CURVE_BUF_SIZE; i++) {
            mCurveLocation[i] = 0;
        }
    }
    
    ~FindCurve() {
        delete mCurveLocation;
    }

    void setThreshold(double threshold) {
        if (threshold == 0)
            threshold = 2;
        mThreshold = threshold;
        mThreshold2 = mThreshold / 2;
    }

    void findCurve() {
        bool start = false;
        bool sign = true;
        const unsigned int N = static_cast<unsigned int>(360 / mThreshold2);
        double * buffer = new double [N];
        double * bufPtr = buffer;
        int index = 0;
        int serial = 0;
        int curveLocationIdx = 0;
        
        for (int i = 0; i < mSize; i++) {
            double angle = mData[i];
            if (!start) {
                if (abs(angle) > mThreshold) {
                    start = true;
                    sign = angle > 0;
                }
            }
            if (start) {
                if (abs(angle) < mThreshold2 || (angle > 0) != sign) {
                    start = false;
                    Curve c;
                    bool ret = c.setData(buffer, index);
                    if (ret) {
                        printf("\nCurve %d\n", serial);
                        c.writeSerial(serial);
                        mCurveLocation[serial] = i - (index - c.mPeak);
                        serial++;
                    }

                    bufPtr = buffer;
                    index = 0;
                } else {
                    *bufPtr++ = angle;
                    index++;
                    if (index >= N) {
                        printf("Run out of buffer\n");
                        return;
                    }
                }
            }
        }

        delete [] buffer;
    }

    void markCurve() {
        Graph<double> g;
        g.setData(mData, mSize);
        g.setThreshold(mThreshold);
        g.setVerticalLines(mCurveLocation);
        g.write("/tmp/pbm/FindCurve.pbm", mSize, 300);
    }
};

int main(int argc, char ** argv) {
    if (argc != 2) {
        puts("Usage:\n\tparser <gps file>\n");
        return -1;
    }

    FILE* fp = fopen(argv[1], "r");
    if (fp == NULL) {
        puts("File open fail\n");
        return -1;
    }

    readFile(fp);

    if (lines == 0) {
        fclose(fp);
        puts("lines is 0\n");
        return -1;
    }

    Vector * vectors = new Vector[lines];

    time_t t0, t1;
    struct tm newyear;
    double seconds;

    time(&t0);  /* get current time; same as: now = time(NULL)  */

    // This can get from the Vincenty formula for 24°57'N 121°33'E
    static float scaleSecondLatitude = 30.887400166;
    static float scaleSecondLongitude = 28.0049474396;
    static float scaleMinuteLatitude = 1853.24878;
    static float scaleMinuteLongitude = 1680.29664;

    struct LatitudeLongitudeSpeed * ptr1 = data, * ptr2 = data + 1;
    for (int i = 0; i < lines - 1; i++) {
		int latDeg = 0, lonDeg = 0;
		if (ptr2->latitude.degree != ptr1->latitude.degree)
			latDeg = (ptr2->latitude.degree - ptr1->latitude.degree) * 3600;
		if (ptr2->longitude.degree != ptr1->longitude.degree)
			lonDeg = (ptr2->longitude.degree - ptr1->longitude.degree) * 3600;

        float distLat = static_cast<float>(ptr2->latitude.minute - ptr1->latitude.minute) * scaleMinuteLatitude +
            (ptr2->latitude.second - ptr1->latitude.second + latDeg) * scaleSecondLatitude;
        float distLon = static_cast<float>(ptr2->longitude.minute - ptr1->longitude.minute) * scaleMinuteLongitude +
            (ptr2->longitude.second - ptr1->longitude.second + lonDeg) * scaleSecondLongitude;

        vectors[i].set(distLon, distLat);
        ptr1++;
        ptr2++;
    }

    double * angles = new double [lines - 2];
    double accx = 0., accy = 0., accl = 0.;
    for (int i = 0; i < lines - 2; i++) {
        Vector& a = vectors[i];
        Vector& b = vectors[i + 1];
        float turn = Vector::crossProduct(a, b);
        float theta = acosf(Vector::cos0(a, b)) * (turn > 0 ? 1 : -1) / 3.141592 * 180;
        if (isnan(theta))  // too small
            theta = 0;
#if 1
        // Distance per meter
        float dist = a.mNorm;
        if (dist < 6.5)
            theta = 0;
        else
            if (abs(theta / dist * 10) < 3)
                theta = 0;
#endif
        angles[i] = theta;

#define DEBUG_LOAD 0
        if (DEBUG_LOAD) {
//            printf("lat:"); data[i].latitude.toString(); printf("\t");
//            printf("lon:"); data[i].longitude.toString(); printf("\t");
            printf("%04d: ", i);
            a.toString();
            printf(",%f,%c\n",
                theta,
                a.mNorm < 3 ? ':' :
                    theta > 3 ? '\\' : theta < -3 ? '/' : '|');
        }
    }
//    Wavelet w(lines - 2);
//    w.setShift(0, 0);
//    w.setShift(1, 1);
//    w.setShift(2, 0);
//    w.setShift(3, 1);
//    w.setData(angles, 0, 6);

    Gaussian g(angles, lines - 2);
    g.filter(8);

    FindCurve f(g.mResult, g.mSize);
    f.findCurve();
    f.markCurve();

    printf("Total %fm\n", accl);

    time(&t1);
    seconds = difftime(t1, t0);
    printf ("%.f seconds for processing %d data.\n", seconds, lines);
    fclose(fp);

    delete [] angles;
    delete [] data;
    delete [] vectors;
    return 0;
}
