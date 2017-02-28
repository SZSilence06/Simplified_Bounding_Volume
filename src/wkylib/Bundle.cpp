#include "Bundle.h"

namespace  WKYLIB {
    Bundle::Bundle()
    {

    }

    void Bundle::checkBeforePut(const std::string &key, KeyValueType valueType)
    {
        auto iterKey = mKeyValueTypes.find(key);
        if(iterKey != mKeyValueTypes.end())
        {
            KeyValueType valueTypeCur = iterKey->second;
            if(valueTypeCur != valueType)
            {
               removeKey(key);
               mKeyValueTypes.insert(std::make_pair(key, valueType));
            }
        }
        else
        {
            mKeyValueTypes.insert(std::make_pair(key, valueType));
        }
    }

    void Bundle::removeKey(const std::string &key)
    {
        auto iterKey = mKeyValueTypes.find(key);
        if(iterKey != mKeyValueTypes.end())
        {
            KeyValueType valueType = iterKey->second;
            switch(valueType)
            {
            case KeyValueType::INT:
                mIntValues.erase(key);
                break;
            case KeyValueType::BOOL:
                mBoolValues.erase(key);
                break;
            case KeyValueType::FLOAT:
                mFloatValues.erase(key);
                break;
            case KeyValueType::DOUBLE:
                mDoubleValues.erase(key);
                break;
            case KeyValueType::STRING:
                mStringValues.erase(key);
                break;
            }
            mKeyValueTypes.erase(iterKey);
        }
    }

    bool Bundle::hasKey(const std::string &key) const
    {
        return mKeyValueTypes.count(key) != 0;
    }

    void Bundle::putInt(const std::string &key, int value)
    {
        checkBeforePut(key, KeyValueType::INT);
        mIntValues.insert(std::make_pair(key, value));
    }

    void Bundle::putBool(const std::string &key, bool value)
    {
        checkBeforePut(key, KeyValueType::BOOL);
        mBoolValues.insert(std::make_pair(key, value));
    }

    void Bundle::putFloat(const std::string &key, float value)
    {
        checkBeforePut(key, KeyValueType::FLOAT);
        mFloatValues.insert(std::make_pair(key, value));
    }

    void Bundle::putDouble(const std::string &key, double value)
    {
        checkBeforePut(key, KeyValueType::DOUBLE);
        mDoubleValues.insert(std::make_pair(key, value));
    }

    void Bundle::putString(const std::string &key, const std::string &value)
    {
        checkBeforePut(key, KeyValueType::STRING);
        mStringValues.insert(std::make_pair(key, value));
    }

    Bundle::KeyValueType Bundle::getKeyValueType(const std::string &key) const
    {
        auto iter = mKeyValueTypes.find(key);
        if(iter == mKeyValueTypes.end())
        {
            throw std::runtime_error("Cannot find KeyValueType when looking for key '" + key + "' in bundle.");
        }
        return iter->second;
    }

    bool Bundle::getInt(const std::string &key, int& value) const
    {
        auto iter = mIntValues.find(key);
        if(iter == mIntValues.end())
        {
            return false;
        }
        value = iter->second;
        return true;
    }

    bool Bundle::getBool(const std::string &key, bool& value) const
    {
        auto iter = mBoolValues.find(key);
        if(iter == mBoolValues.end())
        {
            return false;
        }
        value = iter->second;
        return true;
    }

    bool Bundle::getFloat(const std::string &key, float& value) const
    {
        auto iter = mFloatValues.find(key);
        if(iter == mFloatValues.end())
        {
            return false;
        }
        value = iter->second;
        return true;
    }

    bool Bundle::getDouble(const std::string &key, double& value) const
    {
        auto iter = mDoubleValues.find(key);
        if(iter == mDoubleValues.end())
        {
            return false;
        }
        value = iter->second;
        return true;
    }

    bool Bundle::getString(const std::string &key, std::string& value) const
    {
        auto iter = mStringValues.find(key);
        if(iter == mStringValues.end())
        {
            return false;
        }
        value = iter->second;
        return true;
    }
}
