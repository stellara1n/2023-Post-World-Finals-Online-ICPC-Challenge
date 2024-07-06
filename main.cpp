    #include <bits/stdc++.h>
    #include <stdint.h>
    using namespace std;
    typedef long long ll;

    constexpr int MAXN = 1000005;

    enum SumOp {
        FP16 = 0,
        FP64 = 1,
        FP32 = 2,
        NONE = 3,
    };

    // 实现半精度浮点数
    class Float16{
        static const uint32_t mantissaShift = 42;
        static const uint32_t expShiftMid   = 56;
        static const uint32_t expShiftOut   = 52;
        double dValue_;

    public:
        Float16(double in) : dValue_(in) {
            uint64_t utmp;
            memcpy(&utmp, &dValue_, sizeof utmp);
            //zeroing mantissa bits starting from 11th (this is NOT rounding)
            utmp = utmp >> mantissaShift;
            utmp = utmp << mantissaShift;
            //setting masks for 5-bit exponent extraction out of 11-bit one
            const uint64_t maskExpMid = (63llu << expShiftMid);
            const uint64_t maskExpOut = (15llu << expShiftOut);
            const uint64_t maskExpLead = (1llu << 62);
            const uint64_t maskMantissaD = (1llu << 63) + maskExpLead + maskExpMid + maskExpOut;
            if (utmp & maskExpLead) {// checking leading bit, suspect overflow
                if (utmp & maskExpMid) { //Detected overflow if at least 1 bit is non-zero
                    //Assign Inf with proper sign
                    utmp = utmp | maskExpMid; //setting 1s in the middle 6 bits of of exponent
                    utmp = utmp & maskMantissaD; //zeroing mantissa irrelative of original values to prevent NaN
                    utmp = utmp | maskExpOut; //setting 1s in the last 4 bits of exponent
                }
            } else { //checking small numbers according to exponent range
                if ((utmp & maskExpMid) != maskExpMid) { //Detected underflow if at least 1 bit is 0
                    utmp = 0;
                }
            }
            memcpy(&dValue_, &utmp, sizeof utmp);
        }

        Float16() : dValue_(0) {}

        Float16& operator=(const Float16& rhs) {
            this->dValue_ = rhs.dValue_;
            return *this;
        }

        Float16& operator=(const double& rhs) {
            this->dValue_ = rhs;
            uint64_t utmp;
            memcpy(&utmp, &dValue_, sizeof utmp);
            utmp = utmp >> mantissaShift;
            utmp = utmp << mantissaShift;
            memcpy(&dValue_, &utmp, sizeof utmp);
            return *this;
        }

        friend Float16 operator+(const Float16& lhs, const Float16& rhs) {
            double tmp = lhs.dValue_ + rhs.dValue_;
            return Float16(tmp);
        }

        float convert2Float() { return static_cast<float>(dValue_); }
        double convert2Double() { return dValue_; }
    };

    double calculateFp64(double a, double b) {
        return a+b;
    }

    double calculateFp32(double a, double b) {
        return static_cast<double>(static_cast<float>(a) + static_cast<float>(b));
    }

    double calculateFp16(double a, double b) {
        return (Float16(a) + Float16(b)).convert2Double();
    }

    int d_count = 0, s_count = 0, h_count = 0;

    // 合并
    struct SumNode;
    struct SumNode {
        SumOp op; // 操作
        int idx; // 如果是叶子 记录下标
        double value; // 如果是求和 记录结果
        double exact_value; // 记录准确结果
        int height; // 高度
        SumNode *up; // 上面的节点
        SumNode *left, *right; // 下面的节点
        bool is_local; // 组内的点
        double local_effect; // 组内的影响力
        void print_answer() {
            switch (op) {
                case FP64:
                    d_count++;
                    printf("{d:");
                    left->print_answer();
                    printf(",");
                    right->print_answer();
                    printf("}");
                    break;
                case FP32:
                    s_count++;
                    printf("{s:");
                    left->print_answer();
                    printf(",");
                    right->print_answer();
                    printf("}");
                    break;
                case FP16:
                    h_count++;
                    printf("{h:");
                    left->print_answer();
                    printf(",");
                    right->print_answer();
                    printf("}");
                    break;
                case NONE:
                    printf("%d", idx);
                    break;
                default:
                    printf("error");
                    break;
            }
        }
    };

    struct CompareSumNodePtr {
        bool operator()(const SumNode* a, const SumNode* b) const {
            if (a->height == b->height) {
                return fabs(a->value) > fabs(b->value);
            } else {
                return a->height > b->height;
            }
        }
    };


typedef priority_queue<SumNode*, std::vector<SumNode*>, CompareSumNodePtr> SumNodePriorityQueue;

SumNode* combine_queue(SumNodePriorityQueue& sq, function<SumNode*(SumNode*, SumNode*)> func) {
    if (sq.size() == 0) {
        return nullptr;
    }
    while(sq.size()>1) {
        SumNode* s1 = sq.top();
        sq.pop();
        SumNode* s2 = sq.top();
        sq.pop();
        SumNode* sum = func(s1, s2);
        sq.push(sum);
    }
    SumNode* root = sq.top();
    sq.pop();
    return root;
};


SumNode* merge_global(SumNode* a, SumNode* b) {
    double fp64_value = calculateFp64(a->value, b->value);
    SumNode* merge = new SumNode();
    merge->left = a;
    merge->left->up = merge;
    merge->right = b;
    merge->right->up = merge;
    merge->exact_value = calculateFp64(a->exact_value, b->exact_value);
    merge->op = FP64;
    merge->value = fp64_value;
    merge->height = a->height+1;
    return merge;
}

SumNode* merge_local(SumNode* a, SumNode* b) {
    double fp64_value = calculateFp64(a->value, b->value);
    double fp32_value = calculateFp32(a->value, b->value);
    double fp16_value = calculateFp16(a->value, b->value);
    SumNode* merge = new SumNode();
    merge->left = a;
    merge->left->up = merge;
    merge->right = b;
    merge->right->up = merge;
    merge->exact_value = calculateFp64(a->exact_value, b->exact_value);
    merge->is_local = true;
    merge->height = a->height+1;
    if (!isinf(fp16_value)&&!isnan(fp16_value)) {
        merge->op = FP16;
        merge->value = fp16_value;
    } else if(!isinf(fp32_value)&&!isnan(fp32_value)) {
        merge->op = FP32;
        merge->value = fp32_value;
    } else {
        merge->op = FP64;
        merge->value = fp64_value;
    };
    return merge;
};

// 判题器求和
double kahan_sum(double* arr, int len) {
    double* tmp = new double[len];
    memcpy(tmp, arr, sizeof(double)*len);
    long double trueSum=0, corr=0;
    sort(tmp, tmp+len, [](const double x, const double y) {
        return fabs(x) < fabs(y);
    });
    for (int i=0;i<len;i++) {
        double x = tmp[i];
        volatile long double y = static_cast<long double>(x) - corr;
        volatile long double t = trueSum + y;
        corr = (t - trueSum) - y;
        trueSum = t;
    }
    return (double)trueSum;
}

// 误差公式
double same_mantissa_count(const double calculatedSum, const double expectedSum) {
    if (isinf(calculatedSum) || isinf(expectedSum)) return (double)1.0e20; //worst case
    if (isnan(calculatedSum) || isnan(expectedSum)) return (double)1.0e20; //worst case
    if (calculatedSum == expectedSum) return 1.0e-20; //best case
    return std::max(std::fabs(calculatedSum-expectedSum)/std::max(std::fabs(expectedSum), 1.0e-200), 1.0e-20);
}

// 是否没有误差
bool is_error_tolerable(double se, double sc) {
    return same_mantissa_count(sc, se) < 1.0e-20;
}


int n;
double input[MAXN];
double se = 0.0, sc = 0.0;
SumNode* root;


SumNode* end_nodes[MAXN];
int end_nodes_cnt;

int sorted_idx[MAXN];

SumNode* build_sum_tree_small() {
    for(int i=1;i<=n;i++) {
        sorted_idx[i-1] = i;
    }
    sort(sorted_idx, sorted_idx+n, [](int a, int b) {
        return fabs(input[a]) < fabs(input[b]);
    });
    end_nodes_cnt = 0;
    SumNodePriorityQueue local_queue, body_queue, tail_queue;
    for(int i=0;i<n;i++) {
        int idx = sorted_idx[i];
        SumNode* s = new SumNode();
        s->left = nullptr;
        s->right = nullptr;
        s->op = NONE;
        s->value = input[idx];
        s->exact_value = input[idx];
        s->idx = idx;
        body_queue.push(s);
        end_nodes[end_nodes_cnt++] = s;
    }
    root = combine_queue(body_queue, merge_global);
    return root;
}

SumNode* build_sum_tree(int d) {
    for(int i=1;i<=n;i++) {
        sorted_idx[i-1] = i;
    }
    sort(sorted_idx, sorted_idx+n, [](int a, int b) {
        return fabs(input[a]) < fabs(input[b]);
    });
    int divide_num = n - d;
    
    sort(sorted_idx, sorted_idx + divide_num, [](int a, int b) {
        return a < b;
    });
    end_nodes_cnt = 0;
    SumNodePriorityQueue local_queue, body_queue, tail_queue;
    int m = n/16*16;
    for(int i=0;i<m;i+=16) {
        for(int j=0;j<16;j++) {
            int idx = sorted_idx[i+j];
            SumNode* s = new SumNode();
            s->left = nullptr;
            s->right = nullptr;
            s->op = NONE;
            s->value = input[idx];
            s->exact_value = input[idx];
            s->idx = idx;
            local_queue.push(s);
            end_nodes[end_nodes_cnt++] = s;
        }
        auto local_root = combine_queue(local_queue, merge_local);
        body_queue.push(local_root);
    }
    for(int i=m;i<n;i++) {
        int idx = sorted_idx[i];
        SumNode* s = new SumNode();
        s->left = nullptr;
        s->right = nullptr;
        s->op = NONE;
        s->value = input[idx];
        s->exact_value = input[idx];
        s->idx = idx;
        tail_queue.push(s);
        end_nodes[end_nodes_cnt++] = s;
    }
    // 合并剩余部分
    SumNode *global_root = combine_queue(body_queue, merge_global);
    SumNode *tail_root = combine_queue(tail_queue, merge_local);
    root = nullptr;
    if (global_root == nullptr) {
        root = tail_root;
    } else if (tail_root == nullptr) {
        root = global_root;
    } else {
        root = merge_global(global_root, tail_root);
    }
    return root;
}

int deep = 0;
vector<vector<SumNode*> > deep_vec(MAXN);

void make_deep_info(SumNode* now, int d) {
    deep_vec[d].push_back(now);
    if (now->op==NONE) {
        deep = max(deep, d);
    } else {
        make_deep_info(now->left, d+1);
        make_deep_info(now->right, d+1);
    }
}


int main(int argc, const char * argv[]) {
#ifdef __APPLE__
    freopen("00", "r", stdin);
    freopen("00.ans", "w", stdout);
#endif
    cin >> n;
    for (int i=1; i<=n; i++) {
        cin >> input[i];
    }
    // 准确值
    se = kahan_sum(input+1, n);
    cerr<<setprecision(64)<<"se:"<<se<<endl;
    
    
    // 建树
    for(int i=10000;i<10001;i++) {
        if (n < 20000) {
            root = build_sum_tree_small();
        } else {
            root = build_sum_tree(8000);
        }
        // 实际值
        sc = root->value;
        cerr<<setprecision(64)<<"sc:"<<sc<<endl;
        cerr<<setprecision(64)<<"mantissa-0:"<<same_mantissa_count(sc, se)<<endl;
        for(int i=0;i<end_nodes_cnt;i++) {
            auto node = end_nodes[i];
            auto now = node;
            double tmp = now->value;
            while(now->up != nullptr && now->up->is_local) {
                if (now == now->up->left) {
                    tmp += now->up->right->value;
                } else {
                    tmp += now->up->left->value;
                }
                now = now->up;
            }
            node->local_effect = tmp - now->value;
        }
        
        sort(end_nodes, end_nodes+n, [](SumNode* a, SumNode* b){
            return fabs(se-sc) - fabs(a->local_effect) < fabs(se-sc) - fabs(b->local_effect);
        });
        
        
        for(int i=0;i<end_nodes_cnt;i++) {
            auto node = end_nodes[i];
            if (is_error_tolerable(se, sc)) {
                break;
            }
            if (fabs(sc - se + node->local_effect) < fabs(sc - se)) {
                auto now = node->up;
                while(now!=nullptr) {
                    now->op = FP64;
                    now->value = calculateFp64(now->left->value, now->right->value);
                    now = now->up;
                }
                sc = root->value;
            }
        }
        cerr<<setprecision(64)<<"mantissa-1:"<<same_mantissa_count(sc, se)<<endl;
        
        
        for(int i=1;i<=n;i++) {
            deep_vec[i].resize(0);
        }
        make_deep_info(root, 1);
        for(int i=deep;i>=1;i--) {
            for(int j=0;j<deep_vec[i].size();j++) {
                auto now = deep_vec[i][j];
                if (now->is_local) continue;
                if (now->op == NONE) continue;
                auto old_op = now->op;
                double fp32_value = calculateFp32(now->left->value, now->right->value);
                double fp16_value = calculateFp16(now->left->value, now->right->value);
                if (fabs(sc - se - now->value + fp16_value) <= fabs(sc - se) || is_error_tolerable(se, sc - now->value + fp16_value) ) {
                    now->op = FP16;
                    now->value = fp16_value;
                } else if (fabs(sc - se - now->value + fp32_value) <= fabs(sc - se) || is_error_tolerable(se, sc - now->value + fp32_value) ) {
                    now->op = FP32;
                    now->value = fp32_value;
                }
                if (old_op != now->op) {
                    now = now->up;
                    while(now!=nullptr) {
                        now->value = now->left->value + now->right->value;
                        now = now->up;
                    }
                    sc = root->value;
                }
            }
        }
        cerr<<setprecision(64)<<"mantissa-2:"<<same_mantissa_count(sc, se)<<endl;
    }
    
    root->print_answer();
#ifdef __APPLE__
    cerr<<"count [d:"<<d_count<<", s:"<<s_count<<", h:"<<h_count<<"]"<<endl;
#endif
    
    return 0;
}

