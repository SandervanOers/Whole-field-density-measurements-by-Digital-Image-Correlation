#include "csvload.hpp"
/*--------------------------------------------------------------------------*/
static void removeCharsFromString( std::string &str, char* charsToRemove ) {
   for ( unsigned int i = 0; i < strlen(charsToRemove); ++i ) {
      str.erase( remove(str.begin(), str.end(), charsToRemove[i]), str.end() );
   }
}
/*--------------------------------------------------------------------------*/
class CSVRow
{
    public:
        std::string const& operator[](std::size_t index) const
        {
            return m_data[index];
        }
        std::size_t size() const
        {
            return m_data.size();
        }
        void readNextRow(std::istream& str)
        {
            std::string         line;
            std::getline(str, line);

            std::stringstream   lineStream(line);
            std::string         cell;

            m_data.clear();
            while(std::getline(lineStream, cell, ','))
            {
                m_data.push_back(cell);
            }
            // This checks for a trailing comma with no data after it.
            if (!lineStream && cell.empty())
            {
                // If there was a trailing comma then add an empty element.
                m_data.push_back("");
            }
        }
    private:
        std::vector<std::string> m_data;
};
/*--------------------------------------------------------------------------*/
std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}
/*--------------------------------------------------------------------------*/
class CSVIterator
{
    public:
        typedef std::input_iterator_tag     iterator_category;
        typedef CSVRow                      value_type;
        typedef std::size_t                 difference_type;
        typedef CSVRow*                     pointer;
        typedef CSVRow&                     reference;

        CSVIterator(std::istream& str)  :m_str(str.good()?&str:NULL) { ++(*this); }
        CSVIterator()                   :m_str(NULL) {}

        // Pre Increment
        CSVIterator& operator++()               {if (m_str) { if (!((*m_str) >> m_row)){m_str = NULL;}}return *this;}
        // Post increment
        CSVIterator operator++(int)             {CSVIterator    tmp(*this);++(*this);return tmp;}
        CSVRow const& operator*()   const       {return m_row;}
        CSVRow const* operator->()  const       {return &m_row;}

        bool operator==(CSVIterator const& rhs) {return ((this == &rhs) || ((this->m_str == NULL) && (rhs.m_str == NULL)));}
        bool operator!=(CSVIterator const& rhs) {return !((*this) == rhs);}
    private:
        std::istream*       m_str;
        CSVRow              m_row;
};
/*--------------------------------------------------------------------------*/
extern cv::Mat load_matrix(std::string path, std::string filename, const int &skiplines)
{
  std::ifstream file(path+"/"+filename+".csv");
  std::ifstream file1(path+"/"+filename+".csv");
	int line = 0;
	unsigned int numberofrows, numberofcols;
    for(CSVIterator loop(file); loop != CSVIterator(); ++loop)
    {
  		if (line>=skiplines)
  		{
  			numberofrows = static_cast<unsigned int>((*loop).size());
  		}
		  line++;
    }
	numberofcols = static_cast<unsigned int>(line-skiplines);
	cv::Mat In(numberofcols, numberofrows, CV_64FC1);
	line = 0;
    for(CSVIterator loop(file1); loop != CSVIterator(); ++loop)
    {
		if (line>=skiplines)
		{
			for (unsigned int i = 0; i < (*loop).size(); i++ )
			{
				std::string s = (*loop)[i];
				removeCharsFromString( s, ";" );
				In.at<double>(line-skiplines,i) =  std::stod(s);
			}
		}
		line++;
    }
	return In;
}
/*--------------------------------------------------------------------------*/
