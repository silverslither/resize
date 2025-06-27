import { defineConfig } from "eslint/config";
import globals from "globals";
import path from "node:path";
import { fileURLToPath } from "node:url";
import js from "@eslint/js";
import { FlatCompat } from "@eslint/eslintrc";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const compat = new FlatCompat({
    baseDirectory: __dirname,
    recommendedConfig: js.configs.recommended,
    allConfig: js.configs.all
});

export default defineConfig([{
    extends: compat.extends("eslint:recommended"),
    languageOptions: {
        globals: {
            ...globals.browser,
        },
        ecmaVersion: "latest",
        sourceType: "module",
    },
    rules: {
        indent: ["error", 4, {
            SwitchCase: 1,
        }],
        quotes: ["error", "double"],
        semi: ["error", "always"],
        "prefer-template": ["error"],
        "prefer-const": ["error"],
        "no-trailing-spaces": ["error"],
        "no-multi-spaces": ["error"],
        "object-curly-spacing": ["error", "always"],
        "template-curly-spacing": ["error", "never"],
    },
}]);
