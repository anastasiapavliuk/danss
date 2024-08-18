#pragma once

#ifndef _option_parser_h_included_
#define _option_parser_h_included_

/**
 * Yet another command line options parser
 * Author: Mikhail Veltishchev <dichlofos-mv@yandex.ru>
 **/

#include <string>
#include <iostream>
#include <vector>
#include <deque>
#include <stdlib.h>

namespace NGetOpt {

enum EOptionType {
    ET_STRING,
    ET_INTEGER,
    ET_BOOL,
};

std::string Quote(const std::string& value) {
    return '\'' + value + '\'';
}

struct TOption {
    EOptionType Type;
    char Short;
    std::string Long;
    std::string Description;
    std::vector<std::string> PossibleValues;

    std::string* StrValue;
    int* IntValue;
    bool* BoolValue;

    std::string DefStrValue;
    int DefIntValue;
    bool DefBoolValue;

    TOption()
        : Short(0)
        , StrValue(0)
        , IntValue(0)
        , BoolValue(0)
    {
    }
};

class TOptParser {
public:
    TOptParser(int argc, char** argv) {
        ProgramName = argv[0];
        for (int i = 1; i < argc; ++i) {
            Params.push_back(argv[i]);
        }
    }
    /**
     * Declares an option
     **/
    template<typename T>
    void DeclareOption(char shortOpt, const std::string& longOpt, const std::string& description, T& value, const T& defValue = T(), const std::string& possibleValues = std::string());

    void Parse() {
        // Initialize default values
        for (size_t i = 0; i < Options.size(); ++i) {
            SetDefaultOpt(Options[i]);
        }

        // Process command line arguments
        std::deque<std::string> params;
        for (size_t i = 0; i < Params.size(); ++i) {
            params.push_back(Params[i]);
        };

        while (params.size()) {
            const std::string& arg = params.front();
            if (arg.size() >= 2 && arg[0] == '-' && arg[1] != '-') {
                char shortOpt = arg[1];
                const TOption* opt = Find(shortOpt);
                if (opt) {
                    params.pop_front();
                    ParseOpt(*opt, params);
                    continue;
                }
            }
            if (arg.size() >= 3 && arg.substr(0, 2) == "--") {
                std::string longOpt = arg.substr(2);
                const TOption* opt = Find(longOpt);
                if (opt) {
                    params.pop_front();
                    ParseOpt(*opt, params);
                    continue;
                }
            }
            OnInvalidOption(arg);
        }
    }

private:
    TOption& DeclareOptionImpl(EOptionType type, char shortOpt, const std::string& longOpt, const std::string& description) {
        TOption opt;
        opt.Type = type;
        opt.Short = shortOpt;
        opt.Long = longOpt;
        opt.Description = description;
        Options.push_back(opt);
        return Options.back();
    }

    const TOption* Find(char shortOpt) const {
        for (size_t i = 0; i < Options.size(); ++i) {
            if (Options[i].Short == shortOpt) {
                return &Options[i];
            }
        }
        return 0;
    }

    const TOption* Find(std::string longOpt) const {
        for (size_t i = 0; i < Options.size(); ++i) {
            if (Options[i].Long == longOpt) {
                return &Options[i];
            }
        }
        return 0;
    }

    void ParseOpt(const TOption& opt, std::deque<std::string>& params) {
        switch (opt.Type) {
            case ET_INTEGER: {
                if (!params.size()) {
                    OnMissingValue(opt);
                }
                *opt.IntValue = atoi(params.front().c_str());
                params.pop_front();
                break;
            }
            case ET_STRING: {
                if (!params.size()) {
                    OnMissingValue(opt);
                }
                const std::string& value = params.front();
                if (opt.PossibleValues.size()) {
                    bool found = false;
                    for (size_t i = 0; i < opt.PossibleValues.size(); ++i) {
                        if (value == opt.PossibleValues[i]) {
                            found = true;
                        }
                    }
                    if (!found) {
                        OnInvalidValue(opt);
                    }
                }
                *opt.StrValue = params.front();
                params.pop_front();
                break;
            }
            case ET_BOOL: {
                *opt.BoolValue = true;
                break;
            }
        }
    }

    void SetDefaultOpt(const TOption& opt) {
        switch (opt.Type) {
            case ET_INTEGER: {
                *opt.IntValue = opt.DefIntValue;
                break;
            }
            case ET_STRING: {
                *opt.StrValue = opt.DefStrValue;
                break;
            }
            case ET_BOOL: {
                *opt.BoolValue = opt.DefBoolValue;
                break;
            }
        }
    }

    void OnMissingValue(const TOption& opt) {
        std::cerr << "Option " << Quote("--" + opt.Long)<< " requires a value" << std::endl;
        exit(1);
    }

    void OnInvalidOption(const std::string& arg) {
        std::cerr << "Unknown option " << Quote(arg) << std::endl;
        exit(1);
    }

    void OnInvalidValue(const TOption& opt) {
        std::cerr << "Option " << Quote("--" + opt.Long)<< " can have only these values: ";
        for (size_t i = 0; i < opt.PossibleValues.size(); ++i) {
            if (i) {
                std::cerr << ", ";
            }
            std::cerr << Quote(opt.PossibleValues[i]);
        }
        std::cerr << "." << std::endl;
        exit(1);
    }

private:
    std::string ProgramName;
    std::vector<std::string> Params;
    std::vector<TOption> Options;
};

template<>
void TOptParser::DeclareOption(char shortOpt, const std::string& longOpt, const std::string& description, int& value, const int& defValue, const std::string& possibleValues) {
    TOption& opt = DeclareOptionImpl(ET_INTEGER, shortOpt, longOpt, description);
    opt.DefIntValue = defValue;
    opt.IntValue = &value;
}

template<>
void TOptParser::DeclareOption(char shortOpt, const std::string& longOpt, const std::string& description, std::string& value, const std::string& defValue, const std::string& possibleValues) {
    TOption& opt = DeclareOptionImpl(ET_STRING, shortOpt, longOpt, description);
    std::string current;
    for (size_t i = 0; i < possibleValues.size(); ++i) {
        if (possibleValues[i] == '|') {
            opt.PossibleValues.push_back(current);
            current.clear();
        } else {
            current.push_back(possibleValues[i]);
        }
    }
    if (current.size()) {
        opt.PossibleValues.push_back(current);
    }
    opt.DefStrValue = defValue;
    opt.StrValue = &value;
}

template<>
void TOptParser::DeclareOption(char shortOpt, const std::string& longOpt, const std::string& description, bool& value, const bool& defValue, const std::string& possibleValues) {
    TOption& opt = DeclareOptionImpl(ET_BOOL, shortOpt, longOpt, description);
    opt.BoolValue = &value;
    opt.DefBoolValue = defValue;
}

} // NGetOpt namespace

#endif
