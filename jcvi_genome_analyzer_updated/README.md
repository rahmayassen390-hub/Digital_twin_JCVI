# JCVI Genome Analyzer with AI Agents

تطبيق سطح مكتب احترافي لتحليل الجينوم مع عملاء ذكاء اصطناعي للتحليل المتقدم.

## 📁 هيكل الملفات

```
jcvi_genome_analyzer/
├── __init__.py              # Package initialization
├── main.py                  # Entry point - تشغيل التطبيق
├── data_structures.py       # Data classes (Gene, BlastHit, etc.)
├── neural_networks.py       # AI models (PromoterCNN, TranscriptionLSTM)
├── ai_promoter_agent.py     # AI Agent #1: Promoter Analyzer
├── ai_transcription_agent.py # AI Agent #2: Transcription Dynamics
├── blast_manager.py         # BLAST operations manager
├── genome_manager.py        # Genome data loading and management
├── workers.py               # Background worker threads
├── gui_main.py              # Main GUI window and UI setup
└── gui_methods.py           # GUI event handlers and methods
```

## 📦 الملفات بالتفصيل

### 1. `data_structures.py`
يحتوي على تعريفات البيانات الأساسية:
- `Gene` - بيانات الجين
- `PromoterAnnotation` - بيانات المحفز
- `BlastHit` - نتائج BLAST
- `PromoterPrediction` - تنبؤات المحفز
- `TranscriptionState` - حالة النسخ

### 2. `neural_networks.py`
نماذج الشبكات العصبية:
- `PromoterCNN` - شبكة عصبية التفافية للمحفزات
- `TranscriptionLSTM` - شبكة LSTM للنسخ

### 3. `ai_promoter_agent.py`
AI Agent #1 لتحليل المحفزات:
- تحليل تسلسل upstream
- حساب درجات motif
- تصنيف المحفزات

### 4. `ai_transcription_agent.py`
AI Agent #2 لديناميكيات النسخ:
- استخراج الخصائص الجينية
- التنبؤ بحالة النسخ
- تحديد الاعتمادات

### 5. `blast_manager.py`
إدارة عمليات BLAST:
- البحث عن BLAST
- إنشاء قاعدة بيانات
- تشغيل وتحليل BLAST

### 6. `genome_manager.py`
إدارة بيانات الجينوم:
- تحميل FASTA
- تحميل GFF3
- تحميل ملفات BED

### 7. `workers.py`
Worker threads للعمليات الخلفية:
- `LoadGenomeWorker` - تحميل الجينوم
- `BlastSearchWorker` - بحث BLAST
- `AIAnalysisWorker` - تحليل AI

### 8. `gui_main.py`
واجهة المستخدم الرئيسية:
- تخطيط النوافذ
- إنشاء التبويبات
- عناصر التحكم

### 9. `gui_methods.py`
دوال واجهة المستخدم:
- معالجة الأحداث
- عرض البيانات
- تحريك النسخ

## 🚀 التشغيل

```bash
cd jcvi_genome_analyzer
python main.py
```

## 📋 المتطلبات

```bash
pip install PyQt5 biopython pandas numpy scikit-learn torch openpyxl

# BLAST
conda install -c bioconda blast
# أو
sudo apt-get install ncbi-blast+
```

## ✨ المميزات

1. **تحميل الجينوم المرجعي** (JCVI-syn1.0)
   - ملفات FASTA + GFF3
   - تعليقات المحفزات (BED)

2. **تحميل جينوم الاستعلام** (JCVI-syn3.0)
   - ملفات FASTA
   - تعليقات اختيارية

3. **مقارنة BLAST**
   - تحليل كامل للجينوم
   - ربط الجينات

4. **AI Agent #1: تحليل المحفزات**
   - كشف motif
   - تصنيف المحفزات
   - حساب القوة

5. **AI Agent #2: تحريك النسخ**
   - محاكاة في الوقت الحقيقي
   - تبعيات الجينات
   - ألوان حسب الحالة
