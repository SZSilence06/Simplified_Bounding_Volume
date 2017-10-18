#ifndef WKY_BUNDLE_H
#define WKY_BUNDLE_H

#include <map>
#include <set>
#include <vector>

namespace WKYLIB {
    class Bundle
    {
    public:
        enum KeyValueType
        {
            INT,
            BOOL,
            FLOAT,
            DOUBLE,
            STRING
        };

        template<class ValueType>
        using ValueMap = std::map<std::string, ValueType>;

        using ValueTypeMap = ValueMap<KeyValueType>;
        using IntMap = ValueMap<int>;
        using BoolMap = ValueMap<bool>;
        using FloatMap = ValueMap<float>;
        using DoubleMap = ValueMap<double>;
        using StringMap =  ValueMap<std::string>;

        Bundle();

        void putInt(const std::string& key, int value);
        void putBool(const std::string& key, bool value);
        void putFloat(const std::string& key, float value);
        void putDouble(const std::string& key, double value);
        void putString(const std::string& key, const std::string& value);

        void removeKey(const std::string& key);

        bool hasKey(const std::string& key) const;

        KeyValueType getKeyValueType(const std::string& key) const;
        bool getInt(const std::string& key, int& value) const;
        bool getBool(const std::string& key, bool& value) const;
        bool getFloat(const std::string& key, float& value) const;
        bool getDouble(const std::string& key, double& value) const;
        bool getString(const std::string& key, std::string& value) const;

    private:
        void checkBeforePut(const std::string& key, KeyValueType valueType);

    private:
        ValueTypeMap mKeyValueTypes;

        IntMap mIntValues;
        BoolMap mBoolValues;
        FloatMap mFloatValues;
        DoubleMap mDoubleValues;
        StringMap mStringValues;
    };
}

#endif
